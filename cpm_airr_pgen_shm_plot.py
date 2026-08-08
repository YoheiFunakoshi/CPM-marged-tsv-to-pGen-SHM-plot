#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""CPM exact-UMI AIRR TSV -> pGen, SHM, and pGen-SHM plots.

The primary observation unit is one exact UMI family (one AIRR row).  The
``DUPCOUNT`` annotation embedded in ``sequence_id`` is the number of merged
reads supporting that family and is kept as a separate, optional weight.
Input is an IgBLAST AIRR outfmt 19 TSV, not FASTA or Excel.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import datetime as dt
import io
import json
import math
import os
import re
import sys
import textwrap
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import gaussian_kde

try:
    from openpyxl import Workbook
except Exception:  # pragma: no cover - handled by check_setup / runtime error
    Workbook = None


PGEN_EDGES = [1e-5, 1e-10, 1e-15, 1e-20, 1e-25, 1e-30, 1e-35, 1e-40]
PGEN_LABELS = [
    "1e-5~1e-10",
    "1e-10~1e-15",
    "1e-15~1e-20",
    "1e-20~1e-25",
    "1e-25~1e-30",
    "1e-30~1e-35",
    "1e-35~1e-40",
    "1e-40~",
]
SHM_LABELS = [
    "0~2%",
    "2~4%",
    "4~6%",
    "6~8%",
    "8~10%",
    "10~12%",
    "12~14%",
    "14~16%",
    "16~18%",
    "18~20%",
    "20~%",
]
ROI_DEFINITIONS = [
    {
        "name": "beta1_low_shm_high_pgen",
        "condition": "-16 < log10_pgen < -9 and 0 <= SHM < 2",
        "x_min": -16.0,
        "x_max": -9.0,
        "x_min_inclusive": False,
        "x_max_inclusive": False,
        "y_min": 0.0,
        "y_max": 2.0,
        "y_min_inclusive": True,
        "y_max_inclusive": False,
    },
    {
        "name": "low_shm_core",
        "condition": "-15 <= log10_pgen < -10 and 0 <= SHM < 2",
        "x_min": -15.0,
        "x_max": -10.0,
        "x_min_inclusive": True,
        "x_max_inclusive": False,
        "y_min": 0.0,
        "y_max": 2.0,
        "y_min_inclusive": True,
        "y_max_inclusive": False,
    },
    {
        "name": "low_shm_broad",
        "condition": "-16 <= log10_pgen < -5 and 0 <= SHM < 3",
        "x_min": -16.0,
        "x_max": -5.0,
        "x_min_inclusive": True,
        "x_max_inclusive": False,
        "y_min": 0.0,
        "y_max": 3.0,
        "y_min_inclusive": True,
        "y_max_inclusive": False,
    },
    {
        "name": "main_mid_high_shm_reference",
        "condition": "-21 <= log10_pgen < -10 and 5 <= SHM < 15",
        "x_min": -21.0,
        "x_max": -10.0,
        "x_min_inclusive": True,
        "x_max_inclusive": False,
        "y_min": 5.0,
        "y_max": 15.0,
        "y_min_inclusive": True,
        "y_max_inclusive": False,
    },
]
ROW_LEVEL_COLUMNS = [
    "sequence_id",
    "umi",
    "umi_family_count",
    "supporting_read_count",
    "junction_aa",
    "shm",
    "pgen",
    "log10_pgen",
    "junction",
    "v_identity",
    "locus",
    "productive",
    "v_call",
    "j_call",
    "same_xy_umi_family_count",
    "same_xy_supporting_read_count",
]
AA_ALPHABET = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
TRUTHY = {"T", "TRUE", "1", "Y", "YES"}


LogFn = Callable[[str], None]


@dataclass
class AnalysisConfig:
    input_path: Path
    output_dir: Path
    sample: str
    cache_path: Path
    include_supporting_read_outputs: bool = False
    recalculate_pgen: bool = True
    pgen_workers: int = 6
    min_v_align_len: int = 0
    locus: str = "IGH"
    canonical_mode: str = "rg_reference"
    xlim: tuple[float, float] = (-30.0, -5.0)
    ylim: tuple[float, float] = (0.0, 15.0)
    bw_factor: float = 0.8
    prefix: str | None = None


def log_default(message: str) -> None:
    print(message)


def format_plot_title(value: str, width: int = 68) -> str:
    return "\n".join(textwrap.wrap(str(value), width=width, break_long_words=False))


def safe_sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".igblast.airr.tsv", ".airr.tsv", ".tsv", ".txt", ".zip"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    safe = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", name).strip(" ._")
    return safe[:90] or "sample"


def is_truthy(value: str | None) -> bool:
    return str(value or "").strip().upper() in TRUTHY


def parse_cpm_sequence_id(sequence_id: str | None) -> tuple[str, str, int]:
    """Return base ID, 12-nt UMI, and positive DUPCOUNT from a CPM family ID.

    Collapsed CPM input must preserve exactly one ``UMI=`` and one
    ``DUPCOUNT=`` field in the first whitespace-delimited sequence ID.  A
    malformed or uncollapsed input is rejected rather than silently treated as
    a family count of one.
    """

    value = str(sequence_id or "").strip().split(maxsplit=1)[0]
    if not value:
        raise ValueError("missing sequence_id")
    parts = value.split("|")
    base_id = parts[0].strip()
    umi_values = [part[4:] for part in parts[1:] if part.upper().startswith("UMI=")]
    dup_values = [part[9:] for part in parts[1:] if part.upper().startswith("DUPCOUNT=")]
    if len(umi_values) != 1 or len(dup_values) != 1:
        raise ValueError("sequence_id must contain exactly one UMI and one DUPCOUNT annotation")
    umi = umi_values[0].strip().upper()
    if re.fullmatch(r"[ACGT]{12}", umi) is None:
        raise ValueError("UMI must be exactly 12 A/C/G/T bases")
    if re.fullmatch(r"[0-9]+", dup_values[0].strip()) is None:
        raise ValueError("DUPCOUNT must be a positive integer")
    dupcount = int(dup_values[0])
    if dupcount < 1:
        raise ValueError("DUPCOUNT must be >= 1")
    return base_id, umi, dupcount


def is_acgt(seq: str | None) -> bool:
    seq = str(seq or "").strip()
    return bool(seq) and re.fullmatch(r"[ACGTacgt]+", seq) is not None


def aa_ok(aa: str | None) -> bool:
    aa = str(aa or "").strip().upper()
    if not aa:
        return False
    if "*" in aa or "X" in aa:
        return False
    return all(ch in AA_ALPHABET for ch in aa)


def junction_aa_rejection_reason(aa: str | None, canonical_mode: str) -> str | None:
    value = str(aa or "").strip().upper()
    if not value:
        return "drop_missing_junction_aa"
    if "*" in value or "X" in value or not all(ch in AA_ALPHABET for ch in value):
        return "drop_invalid_junction_aa"
    if canonical_mode == "rg_reference":
        return None
    if canonical_mode != "cpm_conservative":
        raise ValueError(f"Unknown canonical mode: {canonical_mode}")
    if not 5 <= len(value) <= 40:
        return "drop_junction_aa_length_outside_5_40"
    if not value.startswith("C"):
        return "drop_junction_aa_not_cysteine_started"
    if not value.endswith(("F", "W")):
        return "drop_junction_aa_not_fw_ended"
    return None


def to_int(value: str | None, default: int = 0) -> int:
    try:
        text = str(value or "").strip()
        if not text:
            return default
        return int(float(text))
    except Exception:
        return default


def to_float(value: str | None) -> float | None:
    try:
        text = str(value or "").strip()
        if not text:
            return None
        return float(text)
    except Exception:
        return None


def pgen_bin_label(pgen: float) -> str:
    if pgen <= 0:
        return "1e-40~"
    for i in range(len(PGEN_EDGES) - 1):
        hi = PGEN_EDGES[i]
        lo = PGEN_EDGES[i + 1]
        if lo <= pgen < hi:
            return PGEN_LABELS[i]
    if pgen < PGEN_EDGES[-1]:
        return "1e-40~"
    return "1e-5~1e-10"


def shm_bin_label(shm: float) -> str:
    if shm < 0:
        shm = 0.0
    if shm >= 20.0:
        return "20~%"
    start = int(math.floor(shm / 2.0)) * 2
    start = max(0, min(18, start))
    return f"{start}~{start + 2}%"


def shm_from_identity(v_identity: float) -> float:
    if v_identity <= 1.2:
        shm = (1.0 - v_identity) * 100.0
    else:
        shm = 100.0 - v_identity
    return max(0.0, shm)


def ungapped_len(seq: str | None) -> int:
    return len(re.sub(r"[-.\s]", "", str(seq or "")))


def open_airr_tsv(path: Path) -> Iterable[dict[str, str]]:
    """Yield AIRR rows from a plain TSV or a ZIP containing one TSV."""
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            tsv_names = [n for n in zf.namelist() if n.lower().endswith((".tsv", ".txt"))]
            if not tsv_names:
                raise ValueError("ZIP does not contain a TSV/TXT file.")
            if len(tsv_names) > 1:
                raise ValueError("ZIP contains multiple TSV/TXT files. Please unzip and choose one.")
            with zf.open(tsv_names[0], "r") as fb:
                text = io.TextIOWrapper(fb, encoding="utf-8-sig", errors="replace", newline="")
                reader = csv.DictReader(text, delimiter="\t")
                if reader.fieldnames is None:
                    raise ValueError("Input TSV has no header.")
                for row in reader:
                    yield row
        return

    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")
        for row in reader:
            yield row


def get_fieldnames(path: Path) -> list[str]:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            tsv_names = [n for n in zf.namelist() if n.lower().endswith((".tsv", ".txt"))]
            if len(tsv_names) != 1:
                raise ValueError("ZIP must contain exactly one TSV/TXT file.")
            with zf.open(tsv_names[0], "r") as fb:
                header = fb.readline().decode("utf-8-sig", errors="replace").rstrip("\n\r")
        return header.split("\t")
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        return handle.readline().rstrip("\n\r").split("\t")


def v_alignment_len(row: dict[str, str], fields: set[str]) -> int:
    if "v_sequence_alignment" in fields:
        return ungapped_len(row.get("v_sequence_alignment"))
    if "v_alignment_start" in fields and "v_alignment_end" in fields:
        start = to_int(row.get("v_alignment_start"), default=0)
        end = to_int(row.get("v_alignment_end"), default=0)
        if start > 0 and end >= start:
            return end - start + 1
    return 0


def load_pgen_cache(cache_path: Path) -> dict[str, float]:
    cache: dict[str, float] = {}
    if not cache_path.exists():
        return cache
    with cache_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            aa = (row.get("junction_aa") or row.get("AA_JUNCTION") or "").strip().upper()
            pgen = to_float(row.get("pgen"))
            if aa and pgen is not None:
                cache[aa] = pgen
    return cache


_PGEN_MODEL = None


def build_olga_pgen_model():
    try:
        import olga
        import olga.generation_probability as generation_probability
        import olga.load_model as load_model
    except Exception as exc:
        raise RuntimeError("OLGA is not installed. Install with: pip install olga") from exc

    olga_dir = Path(olga.__file__).resolve().parent
    model_dir = olga_dir / "default_models" / "human_B_heavy"
    if not model_dir.exists():
        raise RuntimeError(f"OLGA human_B_heavy model directory not found: {model_dir}")

    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(
        str(model_dir / "model_params.txt"),
        str(model_dir / "V_gene_CDR3_anchors.csv"),
        str(model_dir / "J_gene_CDR3_anchors.csv"),
    )
    gen_model = load_model.GenerativeModelVDJ()
    gen_model.load_and_process_igor_model(str(model_dir / "model_marginals.txt"))
    return generation_probability.GenerationProbabilityVDJ(gen_model, genomic_data)


def compute_one_pgen(aa: str, pgen_model) -> float:
    try:
        return float(pgen_model.compute_aa_CDR3_pgen(aa, print_warnings=False))
    except TypeError:
        try:
            return float(pgen_model.compute_aa_CDR3_pgen(aa))
        except Exception:
            return 0.0
    except Exception:
        return 0.0


def init_pgen_worker() -> None:
    global _PGEN_MODEL
    _PGEN_MODEL = build_olga_pgen_model()


def compute_pgen_worker(aa: str) -> tuple[str, float]:
    global _PGEN_MODEL
    if _PGEN_MODEL is None:
        _PGEN_MODEL = build_olga_pgen_model()
    return aa, compute_one_pgen(aa, _PGEN_MODEL)


def write_fresh_pgen_cache(cache_path: Path, values: dict[str, float]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_name(f"{cache_path.name}.tmp")
    with tmp_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["junction_aa", "pgen"])
        for aa in sorted(values):
            writer.writerow([aa, f"{float(values[aa]):.17g}"])
    os.replace(tmp_path, cache_path)


def compute_pgen_for_aas(
    aas: Iterable[str],
    cache_path: Path,
    log: LogFn,
    recalculate: bool = True,
    workers: int = 6,
) -> dict[str, float]:
    unique_aas = sorted(set(aas))
    workers = max(1, int(workers or 1))
    cache: dict[str, float] = {}
    if recalculate:
        todo = unique_aas
        log(f"pGen cache ignored: recalculating all {len(todo):,} unique junction_aa sequences.")
        log("A fresh pGen cache will be written after successful pGen calculation.")
    else:
        cache = load_pgen_cache(cache_path)
        todo = [aa for aa in unique_aas if aa not in cache]
        if not todo:
            log(f"pGen cache hit: {len(cache):,} entries; no new OLGA calculation.")
            return cache
        log(f"pGen cache enabled: {len(cache):,} cached entries; {len(todo):,} new AA sequences to compute.")

    if not todo:
        return cache

    computed: dict[str, float] = {}
    if workers <= 1 or len(todo) == 1:
        pgen_model = build_olga_pgen_model()
        for i, aa in enumerate(todo, start=1):
            computed[aa] = compute_one_pgen(aa, pgen_model)
            if i % 100 == 0 or i == len(todo):
                log(f"pGen computed {i:,}/{len(todo):,} new AA sequences.")
    else:
        workers = min(workers, len(todo))
        log(f"pGen parallel workers: {workers}")
        executor = None
        try:
            executor = concurrent.futures.ProcessPoolExecutor(max_workers=workers, initializer=init_pgen_worker)
            future_to_aa = {executor.submit(compute_pgen_worker, aa): aa for aa in todo}
            for i, future in enumerate(concurrent.futures.as_completed(future_to_aa), start=1):
                aa = future_to_aa[future]
                try:
                    result_aa, pgen = future.result()
                    computed[result_aa] = float(pgen)
                except Exception:
                    computed[aa] = 0.0
                if i % 100 == 0 or i == len(todo):
                    log(f"pGen computed {i:,}/{len(todo):,} new AA sequences.")
        except Exception as exc:
            remaining = [aa for aa in todo if aa not in computed]
            log(f"Parallel pGen failed ({type(exc).__name__}); falling back to 1 worker for {len(remaining):,} remaining sequences.")
            pgen_model = build_olga_pgen_model()
            for i, aa in enumerate(remaining, start=1):
                computed[aa] = compute_one_pgen(aa, pgen_model)
                if i % 100 == 0 or i == len(remaining):
                    log(f"pGen computed {i:,}/{len(remaining):,} fallback AA sequences.")
        finally:
            if executor is not None:
                try:
                    executor.shutdown(wait=False, cancel_futures=True)
                except TypeError:
                    executor.shutdown(wait=False)

    cache.update(computed)
    if recalculate:
        write_fresh_pgen_cache(cache_path, {aa: cache.get(aa, 0.0) for aa in unique_aas})
        log(f"Saved fresh pGen cache: {cache_path}")
    else:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        need_header = not cache_path.exists()
        with cache_path.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            if need_header:
                writer.writerow(["junction_aa", "pgen"])
            for aa in sorted(computed):
                writer.writerow([aa, f"{computed[aa]:.17g}"])
    return cache


def read_and_aggregate(config: AnalysisConfig, log: LogFn):
    fields = set(get_fieldnames(config.input_path))
    required = {"sequence_id", "productive", "junction", "junction_aa", "v_identity"}
    missing = sorted(required - fields)
    if missing:
        raise ValueError(f"Missing required AIRR column(s): {', '.join(missing)}")

    has_vj = "vj_in_frame" in fields
    has_stop = "stop_codon" in fields
    has_locus = "locus" in fields
    has_vlen = bool({"v_sequence_alignment", "v_alignment_start", "v_alignment_end"} & fields)

    stats = Counter()
    aa_family_counts = Counter()
    aa_supporting_read_counts = Counter()
    j_to_shm: dict[str, list[float]] = defaultdict(list)
    j_to_aa: dict[str, Counter] = defaultdict(Counter)
    j_to_family_count = Counter()
    j_to_supporting_read_count = Counter()
    j_to_vlen: dict[str, list[int]] = defaultdict(list)
    row_records: list[dict[str, object]] = []
    input_umis: set[str] = set()
    kept_umis: set[str] = set()

    for row in open_airr_tsv(config.input_path):
        stats["rows_total"] += 1

        sequence_id = (row.get("sequence_id") or "").strip()
        try:
            _base_id, umi, dupcount = parse_cpm_sequence_id(sequence_id)
        except ValueError as exc:
            raise ValueError(
                f"Invalid CPM exact-UMI sequence_id at AIRR row {stats['rows_total']}: {exc}. "
                "Select the .umi_exact.igblast.airr.tsv generated from collapsed CPM FASTA."
            ) from exc
        input_umis.add(umi)
        stats["input_exact_umi_families"] += 1
        stats["input_supporting_reads"] += dupcount

        locus = (row.get("locus") or "").strip()
        if has_locus and config.locus:
            if locus and locus != config.locus:
                stats["drop_non_IGH_locus"] += 1
                continue
            if not locus:
                stats["kept_empty_locus"] += 1

        if not is_truthy(row.get("productive")):
            stats["drop_nonproductive"] += 1
            continue
        if has_vj and not is_truthy(row.get("vj_in_frame")):
            stats["drop_vj_outframe"] += 1
            continue
        if has_stop and is_truthy(row.get("stop_codon")):
            stats["drop_stopcodon"] += 1
            continue

        junction = (row.get("junction") or "").strip().upper()
        if not is_acgt(junction):
            stats["drop_bad_junction_nt"] += 1
            continue

        aa = (row.get("junction_aa") or "").strip().upper()
        aa_rejection = junction_aa_rejection_reason(aa, config.canonical_mode)
        if aa_rejection:
            stats[aa_rejection] += 1
            continue

        v_identity = to_float(row.get("v_identity"))
        if v_identity is None:
            stats["drop_missing_v_identity"] += 1
            continue

        vlen = v_alignment_len(row, fields) if has_vlen else 0
        if config.min_v_align_len > 0:
            if not has_vlen:
                raise ValueError("min V alignment length was set, but no V alignment length column is available.")
            if vlen < config.min_v_align_len:
                stats["drop_short_v_alignment"] += 1
                continue

        shm = shm_from_identity(v_identity)
        aa_family_counts[aa] += 1
        aa_supporting_read_counts[aa] += dupcount
        j_to_shm[junction].append(shm)
        j_to_aa[junction][aa] += 1
        j_to_family_count[junction] += 1
        j_to_supporting_read_count[junction] += dupcount
        if vlen:
            j_to_vlen[junction].append(vlen)
        kept_umis.add(umi)
        stats["kept_exact_umi_families"] += 1
        stats["kept_supporting_reads"] += dupcount
        row_records.append(
            {
                "sequence_id": sequence_id,
                "umi": umi,
                "umi_family_count": 1,
                "supporting_read_count": dupcount,
                "junction_aa": aa,
                "shm": shm,
                "pgen": 0.0,
                "log10_pgen": math.nan,
                "junction": junction,
                "v_identity": v_identity,
                "locus": locus,
                "productive": (row.get("productive") or "").strip(),
                "v_call": (row.get("v_call") or "").strip(),
                "j_call": (row.get("j_call") or "").strip(),
                "same_xy_umi_family_count": 0,
                "same_xy_supporting_read_count": 0,
            }
        )

    stats["input_unique_umi_sequences"] = len(input_umis)
    stats["kept_unique_umi_sequences"] = len(kept_umis)
    stats["kept_unique_junction_nt"] = len(j_to_family_count)
    stats["kept_unique_junction_aa"] = len(aa_family_counts)
    stats["has_vj_in_frame_column"] = int(has_vj)
    stats["has_stop_codon_column"] = int(has_stop)
    stats["has_locus_column"] = int(has_locus)
    stats["has_sequence_id_umi_dupcount_annotations"] = 1
    stats["has_v_alignment_length_data"] = int(has_vlen)

    if stats["kept_exact_umi_families"] == 0:
        raise ValueError("No exact UMI families remain after productive/canonical/IGH filtering.")
    log(f"Exact UMI families kept after filtering: {stats['kept_exact_umi_families']:,}")
    log(f"DUPCOUNT supporting reads represented by kept families: {stats['kept_supporting_reads']:,}")
    log(f"Junction AA filter mode: {config.canonical_mode}")
    log(f"Unique junction(nt) points: {stats['kept_unique_junction_nt']:,}")
    log(f"Unique junction_aa for pGen: {stats['kept_unique_junction_aa']:,}")
    return (
        aa_family_counts,
        aa_supporting_read_counts,
        j_to_shm,
        j_to_aa,
        j_to_family_count,
        j_to_supporting_read_count,
        j_to_vlen,
        row_records,
        stats,
    )


def write_qc_summary(stats: Counter, out_tsv: Path) -> None:
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key in sorted(stats):
            writer.writerow([key, stats[key]])


def write_pgen_bins(
    aa_family_counts: Counter,
    aa_supporting_read_counts: Counter,
    aa_to_pgen: dict[str, float],
    out_tsv: Path,
) -> tuple[list[float], list[float], list[float]]:
    unique_counts = Counter()
    family_counts = Counter()
    supporting_read_counts = Counter()

    for aa, family_count in aa_family_counts.items():
        label = pgen_bin_label(float(aa_to_pgen.get(aa, 0.0)))
        unique_counts[label] += 1
        family_counts[label] += int(family_count)
        supporting_read_counts[label] += int(aa_supporting_read_counts.get(aa, 0))

    total_unique = sum(unique_counts.values()) or 1
    total_families = sum(family_counts.values()) or 1
    total_supporting_reads = sum(supporting_read_counts.values()) or 1
    frac_unique = [unique_counts.get(label, 0) / total_unique for label in PGEN_LABELS]
    frac_families = [family_counts.get(label, 0) / total_families for label in PGEN_LABELS]
    frac_supporting_reads = [supporting_read_counts.get(label, 0) / total_supporting_reads for label in PGEN_LABELS]

    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "bin",
                "count_unique_junction_aa",
                "fraction_unique_junction_aa",
                "count_exact_umi_families",
                "fraction_exact_umi_families",
                "count_dupcount_supporting_reads",
                "fraction_dupcount_supporting_reads",
            ]
        )
        for label, fu, ff, fs in zip(PGEN_LABELS, frac_unique, frac_families, frac_supporting_reads):
            writer.writerow(
                [
                    label,
                    unique_counts.get(label, 0),
                    f"{fu:.12g}",
                    family_counts.get(label, 0),
                    f"{ff:.12g}",
                    supporting_read_counts.get(label, 0),
                    f"{fs:.12g}",
                ]
            )
    return frac_unique, frac_families, frac_supporting_reads


def plot_barh(labels: list[str], fractions: list[float], ylabel: str, xlabel: str, title: str, out_png: Path) -> None:
    y = np.arange(len(labels))[::-1]
    plt.figure(figsize=(8.4, 5.0))
    plt.barh(y, fractions[::-1], color="#2b7fb8")
    plt.yticks(y, labels[::-1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(format_plot_title(title), fontsize=11)
    xmax = max(fractions) * 1.12 if fractions and max(fractions) > 0 else 1.0
    plt.xlim(0, xmax)
    plt.grid(True, axis="x", alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def write_points(
    junctions: Iterable[str],
    j_to_shm: dict[str, list[float]],
    j_to_aa: dict[str, Counter],
    j_to_family_count: Counter,
    j_to_supporting_read_count: Counter,
    j_to_vlen: dict[str, list[int]],
    aa_to_pgen: dict[str, float],
    out_tsv: Path,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    total_families = sum(j_to_family_count.values()) or 1
    total_supporting_reads = sum(j_to_supporting_read_count.values()) or 1
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "junction_nt",
            "junction_aa",
            "exact_umi_family_count",
            "shm_median",
            "pgen",
            "log10_pgen",
            "v_seq_len_median",
            "aa_candidate_count",
            "exact_umi_family_fraction",
            "dupcount_supporting_read_count",
            "dupcount_supporting_read_fraction",
        ])
        for junction in junctions:
            shm_values = j_to_shm.get(junction, [])
            if not shm_values:
                continue
            aa_counter = j_to_aa.get(junction, Counter())
            representative_aa = aa_counter.most_common(1)[0][0] if aa_counter else ""
            pgen = float(aa_to_pgen.get(representative_aa, 0.0))
            log10_pgen = math.log10(pgen) if pgen > 0 else math.nan
            vlen_values = j_to_vlen.get(junction, [])
            vlen_median = int(np.median(vlen_values)) if vlen_values else 0
            family_count = int(j_to_family_count.get(junction, 0))
            supporting_read_count = int(j_to_supporting_read_count.get(junction, 0))
            row = {
                "junction_nt": junction,
                "junction_aa": representative_aa,
                "exact_umi_family_count": family_count,
                "exact_umi_family_fraction": family_count / total_families,
                "dupcount_supporting_read_count": supporting_read_count,
                "dupcount_supporting_read_fraction": supporting_read_count / total_supporting_reads,
                "shm_median": float(np.median(shm_values)),
                "pgen": pgen,
                "log10_pgen": log10_pgen,
                "v_seq_len_median": vlen_median,
                "aa_candidate_count": len(aa_counter),
            }
            rows.append(row)
            writer.writerow([
                row["junction_nt"],
                row["junction_aa"],
                row["exact_umi_family_count"],
                f"{row['shm_median']:.12g}",
                f"{pgen:.17g}",
                "" if math.isnan(log10_pgen) else f"{log10_pgen:.12g}",
                row["v_seq_len_median"],
                row["aa_candidate_count"],
                f"{row['exact_umi_family_fraction']:.12g}",
                row["dupcount_supporting_read_count"],
                f"{row['dupcount_supporting_read_fraction']:.12g}",
            ])
    return rows


def in_range(value: float, left: float, right: float, left_inclusive: bool, right_inclusive: bool) -> bool:
    left_ok = value >= left if left_inclusive else value > left
    right_ok = value <= right if right_inclusive else value < right
    return left_ok and right_ok


def valid_log10_pgen(value: object) -> bool:
    try:
        numeric = float(value)
    except Exception:
        return False
    return not math.isnan(numeric) and math.isfinite(numeric)


def median_or_blank(values: list[float]) -> str:
    if not values:
        return ""
    return f"{float(np.median(values)):.12g}"


def min_or_blank(values: list[float]) -> str:
    if not values:
        return ""
    return f"{float(min(values)):.12g}"


def max_or_blank(values: list[float]) -> str:
    if not values:
        return ""
    return f"{float(max(values)):.12g}"


def summarize_roi_rows(
    rows: list[dict[str, object]],
    roi: dict[str, object],
    *,
    view: str,
    y_key: str,
    junction_key: str,
    point_weight_key: str,
    family_count_key: str | None,
    supporting_read_count_key: str | None,
    valid_denominator: int,
) -> list[object]:
    selected: list[dict[str, object]] = []
    for row in rows:
        if float(row.get("pgen", 0.0)) <= 0 or not valid_log10_pgen(row.get("log10_pgen")):
            continue
        x = float(row["log10_pgen"])
        y = float(row[y_key])
        if (
            in_range(
                x,
                float(roi["x_min"]),
                float(roi["x_max"]),
                bool(roi["x_min_inclusive"]),
                bool(roi["x_max_inclusive"]),
            )
            and in_range(
                y,
                float(roi["y_min"]),
                float(roi["y_max"]),
                bool(roi["y_min_inclusive"]),
                bool(roi["y_max_inclusive"]),
            )
        ):
            selected.append(row)

    xs = [float(row["log10_pgen"]) for row in selected]
    ys = [float(row[y_key]) for row in selected]
    point_count = len(selected)
    fraction = point_count / valid_denominator if valid_denominator else 0.0
    if family_count_key:
        family_count_sum = sum(float(row.get(family_count_key, 0.0)) for row in selected)
    else:
        family_count_sum = sum(float(row.get(point_weight_key, 1.0)) for row in selected)
    if supporting_read_count_key:
        supporting_read_count_sum = sum(float(row.get(supporting_read_count_key, 0.0)) for row in selected)
    else:
        supporting_read_count_sum = sum(float(row.get(point_weight_key, 1.0)) for row in selected)
    unique_junctions = len({str(row.get(junction_key, "")) for row in selected if str(row.get(junction_key, ""))})
    unique_aas = len({str(row.get("junction_aa", "")) for row in selected if str(row.get("junction_aa", ""))})
    unique_xy = len({(float(row["log10_pgen"]), float(row[y_key])) for row in selected})

    return [
        view,
        roi["name"],
        roi["condition"],
        point_count,
        f"{fraction:.12g}",
        f"{family_count_sum:.12g}",
        f"{supporting_read_count_sum:.12g}",
        unique_junctions,
        unique_aas,
        unique_xy,
        median_or_blank(xs),
        median_or_blank(ys),
        min_or_blank(xs),
        max_or_blank(xs),
        min_or_blank(ys),
        max_or_blank(ys),
    ]


def write_roi_summary(row_records: list[dict[str, object]], beta1_points: list[dict[str, object]], out_tsv: Path) -> None:
    row_valid = [row for row in row_records if float(row.get("pgen", 0.0)) > 0 and valid_log10_pgen(row.get("log10_pgen"))]
    beta_valid = [row for row in beta1_points if float(row.get("pgen", 0.0)) > 0 and valid_log10_pgen(row.get("log10_pgen"))]
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "view",
            "roi_name",
            "condition",
            "point_count",
            "fraction_of_valid_points",
            "exact_umi_family_count_sum",
            "dupcount_supporting_read_count_sum",
            "unique_junction_nt",
            "unique_junction_aa",
            "unique_xy",
            "log10_pgen_median",
            "shm_median",
            "log10_pgen_min_observed",
            "log10_pgen_max_observed",
            "shm_min_observed",
            "shm_max_observed",
        ])
        for roi in ROI_DEFINITIONS:
            writer.writerow(
                summarize_roi_rows(
                    row_valid,
                    roi,
                    view="row_level_exact_umi_family",
                    y_key="shm",
                    junction_key="junction",
                    point_weight_key="umi_family_count",
                    family_count_key="umi_family_count",
                    supporting_read_count_key="supporting_read_count",
                    valid_denominator=len(row_valid),
                )
            )
            writer.writerow(
                summarize_roi_rows(
                    beta_valid,
                    roi,
                    view="beta1_unique_junction_median",
                    y_key="shm_median",
                    junction_key="junction_nt",
                    point_weight_key="exact_umi_family_count",
                    family_count_key="exact_umi_family_count",
                    supporting_read_count_key="dupcount_supporting_read_count",
                    valid_denominator=len(beta_valid),
                )
            )


def add_pgen_to_row_records(row_records: list[dict[str, object]], aa_to_pgen: dict[str, float]) -> None:
    xy_family_counts = Counter()
    xy_supporting_read_counts = Counter()
    for row in row_records:
        aa = str(row["junction_aa"])
        pgen = float(aa_to_pgen.get(aa, 0.0))
        log10_pgen = math.log10(pgen) if pgen > 0 else math.nan
        row["pgen"] = pgen
        row["log10_pgen"] = log10_pgen
        if pgen > 0 and not math.isnan(log10_pgen):
            key = (log10_pgen, float(row["shm"]))
            xy_family_counts[key] += 1
            xy_supporting_read_counts[key] += int(row["supporting_read_count"])

    for row in row_records:
        log10_pgen = float(row["log10_pgen"])
        if float(row["pgen"]) > 0 and not math.isnan(log10_pgen):
            key = (log10_pgen, float(row["shm"]))
            row["same_xy_umi_family_count"] = int(xy_family_counts[key])
            row["same_xy_supporting_read_count"] = int(xy_supporting_read_counts[key])
        else:
            row["same_xy_umi_family_count"] = 0
            row["same_xy_supporting_read_count"] = 0


def clean_excel_value(value: object) -> object:
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def write_row_level_excel(row_records: list[dict[str, object]], out_xlsx: Path) -> None:
    if Workbook is None:
        raise RuntimeError("openpyxl is required to write Excel output. Install it with: pip install openpyxl")
    wb = Workbook(write_only=True)
    ws = wb.create_sheet("accepted_exact_umi_families")
    ws.append(ROW_LEVEL_COLUMNS)
    for row in row_records:
        ws.append([clean_excel_value(row.get(column, "")) for column in ROW_LEVEL_COLUMNS])
    wb.save(out_xlsx)


def write_row_level_tsv(row_records: list[dict[str, object]], out_tsv: Path) -> None:
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(ROW_LEVEL_COLUMNS)
        for row in row_records:
            values = []
            for column in ROW_LEVEL_COLUMNS:
                value = clean_excel_value(row.get(column, ""))
                if value is None:
                    value = ""
                values.append(value)
            writer.writerow(values)


def compact_xy_points(
    rows: list[dict[str, object]],
    config: AnalysisConfig,
    weight_key: str | None = None,
) -> list[dict[str, float]]:
    compacted: dict[tuple[float, float], dict[str, float]] = {}
    for row in rows:
        pgen = float(row["pgen"])
        log10_pgen = float(row["log10_pgen"])
        shm = float(row["shm"])
        if (
            pgen <= 0
            or math.isnan(log10_pgen)
            or not (config.xlim[0] <= log10_pgen <= config.xlim[1])
            or not (config.ylim[0] <= shm <= config.ylim[1])
        ):
            continue
        key = (log10_pgen, shm)
        if key not in compacted:
            compacted[key] = {
                "log10_pgen": log10_pgen,
                "shm": shm,
                "same_xy_count": 0.0,
                "plot_weight": 0.0,
            }
        compacted[key]["same_xy_count"] += 1.0
        if weight_key:
            compacted[key]["plot_weight"] += max(0.0, float(row.get(weight_key, 1.0)))
        else:
            compacted[key]["plot_weight"] += 1.0
    return list(compacted.values())


def write_shm_hist(
    points: list[dict[str, object]],
    out_tsv: Path,
    weight_key: str | None = None,
    shm_key: str = "shm_median",
    unit_name: str = "exact_umi_families",
) -> list[float]:
    counts = Counter()
    for row in points:
        weight = float(row.get(weight_key, 1.0)) if weight_key else 1.0
        counts[shm_bin_label(float(row[shm_key]))] += weight
    total = sum(counts.values()) or 1
    fractions = [counts.get(label, 0) / total for label in SHM_LABELS]
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["bin", f"count_{unit_name}", f"fraction_{unit_name}"])
        for label, fraction in zip(SHM_LABELS, fractions):
            count_value = counts.get(label, 0)
            if float(count_value).is_integer():
                count_value = int(count_value)
            writer.writerow([label, count_value, f"{fraction:.12g}"])
    return fractions


def plot_kde(
    points: list[dict[str, object]],
    config: AnalysisConfig,
    out_png: Path,
    log: LogFn,
    weight_key: str | None = None,
    title_suffix: str | None = None,
) -> int:
    filtered = [
        row
        for row in points
        if float(row["pgen"]) > 0
        and not math.isnan(float(row["log10_pgen"]))
        and config.xlim[0] <= float(row["log10_pgen"]) <= config.xlim[1]
        and config.ylim[0] <= float(row["shm_median"]) <= config.ylim[1]
    ]
    if len(filtered) < 5:
        log("Too few valid points for KDE; writing scatter fallback instead.")
        plot_scatter(filtered, config, out_png)
        return len(filtered)

    xs = np.array([float(row["log10_pgen"]) for row in filtered], dtype=float)
    ys = np.array([float(row["shm_median"]) for row in filtered], dtype=float)
    weights = None
    if weight_key:
        weights = np.array([max(0.0, float(row.get(weight_key, 0.0))) for row in filtered], dtype=float)
        if weights.sum() <= 0:
            weights = None
    try:
        kde_args = {
            "bw_method": lambda obj: obj.scotts_factor() * config.bw_factor,
        }
        if weights is not None:
            kde_args["weights"] = weights
        kde = gaussian_kde(np.vstack([xs, ys]), **kde_args)
        xi = np.linspace(config.xlim[0], config.xlim[1], 250)
        yi = np.linspace(config.ylim[0], config.ylim[1], 250)
        x_grid, y_grid = np.meshgrid(xi, yi)
        density = kde(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)
    except TypeError as exc:
        if weights is not None:
            log(f"Weighted KDE is not supported by this SciPy ({exc}); writing unweighted KDE instead.")
            return plot_kde(points, config, out_png, log, weight_key=None, title_suffix=title_suffix)
        log(f"KDE failed ({exc}); writing scatter fallback instead.")
        plot_scatter(filtered, config, out_png)
        return len(filtered)
    except Exception as exc:
        log(f"KDE failed ({exc}); writing scatter fallback instead.")
        plot_scatter(filtered, config, out_png)
        return len(filtered)

    plt.figure(figsize=(6.2, 6.0))
    plt.contourf(x_grid, y_grid, density, levels=12, cmap="YlOrRd")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    suffix = title_suffix or ("weighted KDE" if weight_key else "KDE")
    plt.title(format_plot_title(f"{config.sample} pGen-SHM {suffix}"), fontsize=11)
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return len(filtered)


def plot_weighted_scatter(points: list[dict[str, object]], config: AnalysisConfig, out_png: Path) -> int:
    filtered = [
        row
        for row in points
        if float(row["pgen"]) > 0
        and not math.isnan(float(row["log10_pgen"]))
        and config.xlim[0] <= float(row["log10_pgen"]) <= config.xlim[1]
        and config.ylim[0] <= float(row["shm_median"]) <= config.ylim[1]
    ]
    xs = np.array([float(row["log10_pgen"]) for row in filtered], dtype=float)
    ys = np.array([float(row["shm_median"]) for row in filtered], dtype=float)
    weights = np.array([max(1.0, float(row.get("dupcount_supporting_read_count", 1))) for row in filtered], dtype=float)
    if len(filtered) == 0:
        sizes = np.array([], dtype=float)
    else:
        log_weights = np.log10(weights + 1.0)
        if float(log_weights.max()) > float(log_weights.min()):
            sizes = 12.0 + 80.0 * (log_weights - log_weights.min()) / (log_weights.max() - log_weights.min())
        else:
            sizes = np.full_like(log_weights, 28.0)

    plt.figure(figsize=(6.2, 6.0))
    plt.scatter(xs, ys, s=sizes, alpha=0.38, color="#2b7fb8", edgecolors="none")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(
        format_plot_title(f"{config.sample} pGen-SHM scatter weighted by DUPCOUNT supporting reads"),
        fontsize=11,
    )
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return len(filtered)


def plot_xy_kde(
    xy_points: list[dict[str, float]],
    config: AnalysisConfig,
    out_png: Path,
    log: LogFn,
    title_suffix: str,
    log_density: bool = False,
) -> int:
    if len(xy_points) < 5:
        log("Too few valid row-level points for KDE; writing scatter fallback instead.")
        plot_xy_scatter(xy_points, config, out_png, title_suffix)
        return int(sum(point.get("plot_weight", 1.0) for point in xy_points))

    xs = np.array([float(point["log10_pgen"]) for point in xy_points], dtype=float)
    ys = np.array([float(point["shm"]) for point in xy_points], dtype=float)
    weights = np.array([max(0.0, float(point.get("plot_weight", 1.0))) for point in xy_points], dtype=float)
    if weights.sum() <= 0:
        weights = None
    try:
        kde_args = {
            "bw_method": lambda obj: obj.scotts_factor() * config.bw_factor,
        }
        if weights is not None:
            kde_args["weights"] = weights
        kde = gaussian_kde(np.vstack([xs, ys]), **kde_args)
        xi = np.linspace(config.xlim[0], config.xlim[1], 250)
        yi = np.linspace(config.ylim[0], config.ylim[1], 250)
        x_grid, y_grid = np.meshgrid(xi, yi)
        density = kde(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)
    except TypeError as exc:
        log(f"Weighted KDE is not supported by this SciPy ({exc}); writing scatter fallback instead.")
        plot_xy_scatter(xy_points, config, out_png, title_suffix)
        return int(sum(point.get("plot_weight", 1.0) for point in xy_points))
    except Exception as exc:
        log(f"KDE failed ({exc}); writing scatter fallback instead.")
        plot_xy_scatter(xy_points, config, out_png, title_suffix)
        return int(sum(point.get("plot_weight", 1.0) for point in xy_points))

    plt.figure(figsize=(6.2, 6.0))
    if log_density:
        positive_density = density[density > 0]
        if positive_density.size and float(density.max()) > 0:
            density_max = float(density.max())
            density_min = max(float(np.percentile(positive_density, 1)), density_max * 0.005)
            if density_min < density_max:
                levels = np.geomspace(density_min, density_max, 12)
                plt.contourf(
                    x_grid,
                    y_grid,
                    density,
                    levels=levels,
                    cmap="YlOrRd",
                    norm=LogNorm(vmin=density_min, vmax=density_max),
                    extend="min",
                )
            else:
                plt.contourf(x_grid, y_grid, density, levels=12, cmap="YlOrRd")
        else:
            plt.contourf(x_grid, y_grid, density, levels=12, cmap="YlOrRd")
    else:
        plt.contourf(x_grid, y_grid, density, levels=12, cmap="YlOrRd")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(format_plot_title(f"{config.sample} pGen-SHM {title_suffix}"), fontsize=11)
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return int(sum(point.get("plot_weight", 1.0) for point in xy_points))


def plot_xy_scatter(
    xy_points: list[dict[str, float]],
    config: AnalysisConfig,
    out_png: Path,
    title_suffix: str,
) -> int:
    xs = np.array([float(point["log10_pgen"]) for point in xy_points], dtype=float)
    ys = np.array([float(point["shm"]) for point in xy_points], dtype=float)
    weights = np.array([max(1.0, float(point.get("plot_weight", 1.0))) for point in xy_points], dtype=float)
    if len(xy_points) == 0:
        sizes = np.array([], dtype=float)
    else:
        log_weights = np.log10(weights + 1.0)
        if float(log_weights.max()) > float(log_weights.min()):
            sizes = 10.0 + 75.0 * (log_weights - log_weights.min()) / (log_weights.max() - log_weights.min())
        else:
            sizes = np.full_like(log_weights, 24.0)

    plt.figure(figsize=(6.2, 6.0))
    plt.scatter(xs, ys, s=sizes, alpha=0.35, color="#2b7fb8", edgecolors="none")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(format_plot_title(f"{config.sample} pGen-SHM {title_suffix}"), fontsize=11)
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return int(sum(point.get("plot_weight", 1.0) for point in xy_points))


def plot_scatter(points: list[dict[str, object]], config: AnalysisConfig, out_png: Path) -> None:
    xs = [float(row["log10_pgen"]) for row in points if float(row["pgen"]) > 0]
    ys = [float(row["shm_median"]) for row in points if float(row["pgen"]) > 0]
    plt.figure(figsize=(6.2, 6.0))
    plt.scatter(xs, ys, s=10, alpha=0.4, color="#2b7fb8", edgecolors="none")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(format_plot_title(f"{config.sample} pGen-SHM"), fontsize=11)
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def write_run_log(lines: list[str], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")


def write_run_conditions(config: AnalysisConfig, stats: Counter, out_path: Path) -> None:
    payload = {
        "application": "CPM exact-UMI AIRR TSV to pGen-SHM plot",
        "created_at": dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds"),
        "input_airr_tsv": str(config.input_path),
        "output_directory": str(config.output_dir),
        "sample": config.sample,
        "primary_observation_unit": "exact_umi_family",
        "primary_unit_definition": "one AIRR row with a valid 12-nt UMI and positive DUPCOUNT",
        "secondary_weight_unit": "DUPCOUNT supporting merged reads",
        "supporting_read_weighted_plots_created": bool(config.include_supporting_read_outputs),
        "filters": {
            "locus": config.locus,
            "productive": True,
            "vj_in_frame_if_present": True,
            "stop_codon_if_present": False,
            "junction_nt": "A/C/G/T only",
            "junction_aa": "non-empty; no stop codon or X",
            "junction_aa_filter_mode": config.canonical_mode,
            "min_v_alignment_length": int(config.min_v_align_len),
        },
        "shm_definition": "100 - v_identity_percent; fractional v_identity is converted to percent",
        "pgen_definition": "OLGA human_B_heavy compute_aa_CDR3_pgen(junction_aa)",
        "pgen_cache": str(config.cache_path),
        "pgen_cache_recalculated": bool(config.recalculate_pgen),
        "pgen_workers": int(config.pgen_workers),
        "plot": {
            "xlim_log10_pgen": list(config.xlim),
            "ylim_shm_percent": list(config.ylim),
            "bandwidth_factor": config.bw_factor,
        },
        "qc": dict(sorted(stats.items())),
    }
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def run_analysis(config: AnalysisConfig, log: LogFn = log_default) -> dict[str, Path]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.prefix or safe_sample_name(Path(config.sample))
    outputs = {
        "qc_summary": config.output_dir / f"{prefix}_qc_summary.tsv",
        "pgen_bins": config.output_dir / f"{prefix}_pgen_bins.tsv",
        "pgen_bins_unique_png": config.output_dir / f"{prefix}_pgen_bins_unique_junction_aa.png",
        "pgen_bins_umi_family_png": config.output_dir / f"{prefix}_pgen_bins_exact_umi_family.png",
        "shm_hist_umi_family": config.output_dir / f"{prefix}_shm_hist_exact_umi_family.tsv",
        "shm_hist_umi_family_png": config.output_dir / f"{prefix}_shm_hist_exact_umi_family.png",
        "rows_xlsx": config.output_dir / f"{prefix}_pgen_shm_rows.xlsx",
        "points": config.output_dir / f"{prefix}_pgen_shm_points.tsv",
        "beta1_points": config.output_dir / f"{prefix}_pgen_shm_beta1_unique_junction_points.tsv",
        "roi_summary": config.output_dir / f"{prefix}_pgen_shm_roi_summary.tsv",
        "kde_umi_family_png": config.output_dir / f"{prefix}_pgen_shm_kde_exact_umi_family.png",
        "kde_umi_family_log_density_png": config.output_dir / f"{prefix}_pgen_shm_kde_exact_umi_family_log_density.png",
        "beta1_kde_png": config.output_dir / f"{prefix}_pgen_shm_kde_beta1_unique_junction_unweighted.png",
        "beta1_kde_umi_family_png": config.output_dir / f"{prefix}_pgen_shm_kde_beta1_exact_umi_family.png",
        "scatter_umi_family_png": config.output_dir / f"{prefix}_pgen_shm_scatter_exact_umi_family.png",
        "run_log": config.output_dir / f"{prefix}_run_log.txt",
        "run_conditions": config.output_dir / f"{prefix}_run_conditions.json",
    }
    if config.include_supporting_read_outputs:
        outputs.update(
            {
                "pgen_bins_supporting_reads_png": config.output_dir / f"{prefix}_pgen_bins_dupcount_supporting_reads.png",
                "shm_hist_supporting_reads": config.output_dir / f"{prefix}_shm_hist_dupcount_supporting_reads.tsv",
                "shm_hist_supporting_reads_png": config.output_dir / f"{prefix}_shm_hist_dupcount_supporting_reads.png",
                "kde_supporting_reads_png": config.output_dir / f"{prefix}_pgen_shm_kde_dupcount_supporting_reads.png",
                "kde_supporting_reads_log_density_png": config.output_dir / f"{prefix}_pgen_shm_kde_dupcount_supporting_reads_log_density.png",
                "scatter_supporting_reads_png": config.output_dir / f"{prefix}_pgen_shm_scatter_dupcount_supporting_reads.png",
                "beta1_kde_supporting_reads_png": config.output_dir / f"{prefix}_pgen_shm_kde_beta1_dupcount_supporting_reads.png",
            }
        )
    run_log_lines: list[str] = []

    def log_both(message: str) -> None:
        run_log_lines.append(message)
        log(message)

    log_both(f"Input AIRR TSV: {config.input_path}")
    log_both(f"Output folder: {config.output_dir}")
    log_both(f"Sample: {config.sample}")
    log_both(f"pGen cache: {config.cache_path}")
    log_both(
        f"pGen cache policy: {'recalculate all and refresh cache' if config.recalculate_pgen else 'use existing cache when available'}"
    )
    log_both(f"pGen workers: {max(1, int(config.pgen_workers))}")
    log_both("Primary observation unit: exact UMI family (one AIRR row).")
    log_both("Secondary unit: DUPCOUNT supporting merged reads; never pooled with family counts.")
    log_both(
        f"DUPCOUNT-weighted plot output: {'enabled' if config.include_supporting_read_outputs else 'disabled'}"
    )
    log_both(
        f"Locus policy: keep {config.locus}; empty locus is retained and counted in QC, matching the RG reference."
    )
    log_both(f"Junction AA filter mode: {config.canonical_mode}")

    (
        aa_family_counts,
        aa_supporting_read_counts,
        j_to_shm,
        j_to_aa,
        j_to_family_count,
        j_to_supporting_read_count,
        j_to_vlen,
        row_records,
        stats,
    ) = read_and_aggregate(config, log_both)
    write_qc_summary(stats, outputs["qc_summary"])
    log_both(f"Saved QC summary: {outputs['qc_summary']}")

    aa_all = sorted(aa_family_counts.keys())
    aa_to_pgen = compute_pgen_for_aas(
        aa_all,
        config.cache_path,
        log_both,
        recalculate=config.recalculate_pgen,
        workers=config.pgen_workers,
    )
    add_pgen_to_row_records(row_records, aa_to_pgen)
    beta1_points = write_points(
        sorted(j_to_family_count.keys()),
        j_to_shm,
        j_to_aa,
        j_to_family_count,
        j_to_supporting_read_count,
        j_to_vlen,
        aa_to_pgen,
        outputs["beta1_points"],
    )
    log_both(f"Saved beta1-compatible unique junction points: {outputs['beta1_points']}")

    frac_unique, frac_families, frac_supporting_reads = write_pgen_bins(
        aa_family_counts,
        aa_supporting_read_counts,
        aa_to_pgen,
        outputs["pgen_bins"],
    )
    plot_barh(
        PGEN_LABELS,
        frac_unique,
        "pGen",
        "Fraction of unique junction AA",
        f"{config.sample} pGen bins (unique junction AA)",
        outputs["pgen_bins_unique_png"],
    )
    plot_barh(
        PGEN_LABELS,
        frac_families,
        "pGen",
        "Fraction of exact UMI families",
        f"{config.sample} pGen bins (exact UMI families)",
        outputs["pgen_bins_umi_family_png"],
    )
    if config.include_supporting_read_outputs:
        plot_barh(
            PGEN_LABELS,
            frac_supporting_reads,
            "pGen",
            "Fraction of DUPCOUNT supporting reads",
            f"{config.sample} pGen bins (DUPCOUNT supporting reads)",
            outputs["pgen_bins_supporting_reads_png"],
        )
    log_both(f"Saved pGen bins and primary plots: {outputs['pgen_bins']}")

    write_row_level_excel(row_records, outputs["rows_xlsx"])
    write_row_level_tsv(row_records, outputs["points"])
    zero_points = sum(1 for row in row_records if float(row["pgen"]) <= 0)
    duplicate_xy_rows = sum(
        1 for row in row_records if int(row.get("same_xy_umi_family_count", 0)) > 1
    )
    unique_xy_points = len(
        {
            (float(row["log10_pgen"]), float(row["shm"]))
            for row in row_records
            if float(row["pgen"]) > 0
        }
    )
    log_both(f"Saved family-level pGen-SHM Excel: {outputs['rows_xlsx']}")
    log_both(f"Saved family-level pGen-SHM points TSV: {outputs['points']}")
    log_both(f"pGen=0 families retained in tables and excluded from plots: {zero_points:,}")
    log_both(
        f"Family-level rows: {len(row_records):,}; unique x-y coordinates: {unique_xy_points:,}; "
        f"families sharing an x-y coordinate: {duplicate_xy_rows:,}"
    )

    write_roi_summary(row_records, beta1_points, outputs["roi_summary"])
    log_both(f"Saved pGen-SHM ROI summary: {outputs['roi_summary']}")

    shm_family_frac = write_shm_hist(
        row_records,
        outputs["shm_hist_umi_family"],
        weight_key="umi_family_count",
        shm_key="shm",
        unit_name="exact_umi_families",
    )
    plot_barh(
        SHM_LABELS,
        shm_family_frac,
        "%Mutation",
        "Fraction of exact UMI families",
        f"{config.sample} SHM histogram (exact UMI families)",
        outputs["shm_hist_umi_family_png"],
    )
    log_both(
        f"Saved exact-UMI-family SHM histogram: {outputs['shm_hist_umi_family']} ; "
        f"{outputs['shm_hist_umi_family_png']}"
    )

    family_xy_points = compact_xy_points(row_records, config, weight_key="umi_family_count")
    family_scatter_weight = plot_xy_scatter(
        family_xy_points,
        config,
        outputs["scatter_umi_family_png"],
        "exact UMI family scatter",
    )
    family_kde_weight = plot_xy_kde(
        family_xy_points,
        config,
        outputs["kde_umi_family_png"],
        log_both,
        "exact UMI family KDE",
    )
    family_log_weight = plot_xy_kde(
        family_xy_points,
        config,
        outputs["kde_umi_family_log_density_png"],
        log_both,
        "exact UMI family log-density KDE",
        log_density=True,
    )
    log_both(
        f"Primary plots use exact UMI families: scatter {family_scatter_weight:,}; "
        f"KDE {family_kde_weight:,}; log-density KDE {family_log_weight:,}."
    )

    beta1_kde_points = plot_kde(
        beta1_points,
        config,
        outputs["beta1_kde_png"],
        log_both,
        weight_key=None,
        title_suffix="beta1 unique-junction unweighted KDE",
    )
    beta1_family_kde_points = plot_kde(
        beta1_points,
        config,
        outputs["beta1_kde_umi_family_png"],
        log_both,
        weight_key="exact_umi_family_count",
        title_suffix="beta1 unique-junction exact-UMI-family KDE",
    )
    log_both(
        f"Beta1 views: unweighted unique-junction points {beta1_kde_points:,}; "
        f"family-weighted coordinates {beta1_family_kde_points:,}."
    )

    if config.include_supporting_read_outputs:
        shm_support_frac = write_shm_hist(
            row_records,
            outputs["shm_hist_supporting_reads"],
            weight_key="supporting_read_count",
            shm_key="shm",
            unit_name="dupcount_supporting_reads",
        )
        plot_barh(
            SHM_LABELS,
            shm_support_frac,
            "%Mutation",
            "Fraction of DUPCOUNT supporting reads",
            f"{config.sample} SHM histogram (DUPCOUNT supporting reads)",
            outputs["shm_hist_supporting_reads_png"],
        )
        supporting_xy_points = compact_xy_points(
            row_records, config, weight_key="supporting_read_count"
        )
        support_kde_weight = plot_xy_kde(
            supporting_xy_points,
            config,
            outputs["kde_supporting_reads_png"],
            log_both,
            "DUPCOUNT supporting-read KDE",
        )
        support_log_weight = plot_xy_kde(
            supporting_xy_points,
            config,
            outputs["kde_supporting_reads_log_density_png"],
            log_both,
            "DUPCOUNT supporting-read log-density KDE",
            log_density=True,
        )
        support_scatter_weight = plot_xy_scatter(
            supporting_xy_points,
            config,
            outputs["scatter_supporting_reads_png"],
            "DUPCOUNT supporting-read scatter",
        )
        beta1_support_weight = plot_kde(
            beta1_points,
            config,
            outputs["beta1_kde_supporting_reads_png"],
            log_both,
            weight_key="dupcount_supporting_read_count",
            title_suffix="beta1 unique-junction DUPCOUNT-support KDE",
        )
        log_both(
            "Secondary plots use DUPCOUNT supporting reads: "
            f"scatter {support_scatter_weight:,}; KDE {support_kde_weight:,}; "
            f"log-density KDE {support_log_weight:,}; beta1 coordinates {beta1_support_weight:,}."
        )

    write_run_conditions(config, stats, outputs["run_conditions"])
    log_both(f"Saved run conditions JSON: {outputs['run_conditions']}")
    write_run_log(run_log_lines, outputs["run_log"])
    log(f"Saved run log: {outputs['run_log']}")
    return outputs


def parse_range(value: str, label: str) -> tuple[float, float]:
    try:
        left, right = [float(x.strip()) for x in value.split(",", 1)]
    except Exception as exc:
        raise argparse.ArgumentTypeError(f"{label} must be two comma-separated numbers, e.g. -30,-5") from exc
    if right <= left:
        raise argparse.ArgumentTypeError(f"{label}: max must be greater than min.")
    return left, right


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CPM exact-UMI AIRR TSV -> pGen / SHM / pGen-SHM plots")
    parser.add_argument("--input", required=True, help="CPM .umi_exact IgBLAST AIRR TSV, TXT, or ZIP containing one TSV.")
    parser.add_argument("--outdir", default="", help="Output folder. Default: same folder as input.")
    parser.add_argument("--sample", default="", help="Sample name for plot titles and output prefix.")
    parser.add_argument("--pgen-cache", default="", help="pGen cache TSV. Default: <outdir>/pgen_cache.tsv")
    parser.add_argument("--use-pgen-cache", action="store_true", help="Use existing pGen cache when available. Default: recalculate all pGen.")
    parser.add_argument("--pgen-workers", type=int, default=6, help="Parallel OLGA pGen worker processes. Default: 6.")
    parser.add_argument(
        "--include-supporting-read-outputs",
        action="store_true",
        help="Also create secondary plots weighted by sequence_id DUPCOUNT supporting reads.",
    )
    parser.add_argument("--min-v-align-len", type=int, default=0, help="Optional V alignment length filter. Default: 0 disabled.")
    parser.add_argument("--locus", default="IGH", help="Expected locus. Default: IGH. Empty disables locus filtering.")
    parser.add_argument(
        "--canonical-mode",
        choices=("rg_reference", "cpm_conservative"),
        default="rg_reference",
        help=(
            "junction_aa filter: rg_reference keeps letters without X/stop; "
            "cpm_conservative additionally requires C...F/W and length 5-40."
        ),
    )
    parser.add_argument("--xlim", default="-30,-5", help="KDE x-axis log10 pGen range, e.g. -30,-5.")
    parser.add_argument("--ylim", default="0,15", help="KDE y-axis SHM range, e.g. 0,15.")
    parser.add_argument("--bw-factor", type=float, default=0.8, help="KDE bandwidth multiplier. Default: 0.8.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    input_path = Path(args.input).resolve()
    if not input_path.exists():
        parser.error(f"input not found: {input_path}")
    output_dir = Path(args.outdir).resolve() if args.outdir else input_path.parent
    sample = args.sample.strip() or safe_sample_name(input_path)
    cache_path = Path(args.pgen_cache).resolve() if args.pgen_cache else output_dir / "pgen_cache.tsv"
    config = AnalysisConfig(
        input_path=input_path,
        output_dir=output_dir,
        sample=sample,
        cache_path=cache_path,
        include_supporting_read_outputs=bool(args.include_supporting_read_outputs),
        recalculate_pgen=not bool(args.use_pgen_cache),
        pgen_workers=max(1, int(args.pgen_workers)),
        min_v_align_len=max(0, int(args.min_v_align_len)),
        locus=args.locus.strip(),
        canonical_mode=args.canonical_mode,
        xlim=parse_range(args.xlim, "xlim"),
        ylim=parse_range(args.ylim, "ylim"),
        bw_factor=float(args.bw_factor),
        prefix=safe_sample_name(Path(sample)),
    )
    run_analysis(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
