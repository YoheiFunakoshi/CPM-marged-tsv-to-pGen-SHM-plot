from __future__ import annotations

import csv
import unittest
from pathlib import Path

from cpm_airr_pgen_shm_plot import (
    AnalysisConfig,
    junction_aa_rejection_reason,
    parse_cpm_sequence_id,
    read_and_aggregate,
    shm_from_identity,
)


FIELDS = [
    "sequence_id",
    "productive",
    "vj_in_frame",
    "stop_codon",
    "locus",
    "junction",
    "junction_aa",
    "v_identity",
    "v_sequence_alignment",
    "v_call",
    "j_call",
]


def row(sequence_id: str, *, productive: str = "T", dup_junction: str = "TGTGCCAGG") -> dict[str, str]:
    return {
        "sequence_id": sequence_id,
        "productive": productive,
        "vj_in_frame": "T",
        "stop_codon": "F",
        "locus": "IGH",
        "junction": dup_junction,
        "junction_aa": "CAR",
        "v_identity": "98.5",
        "v_sequence_alignment": "ACGTACGT",
        "v_call": "IGHV1-2*02",
        "j_call": "IGHJ4*02",
    }


class CpmUnitTests(unittest.TestCase):
    def test_parse_valid_sequence_id(self) -> None:
        base, umi, dupcount = parse_cpm_sequence_id(
            "M03847:1:FC:1:1101:10:20|UMI=ACGTACGTACGT|DUPCOUNT=7"
        )
        self.assertEqual(base, "M03847:1:FC:1:1101:10:20")
        self.assertEqual(umi, "ACGTACGTACGT")
        self.assertEqual(dupcount, 7)

    def test_parse_rejects_missing_or_malformed_annotations(self) -> None:
        invalid = [
            "read1",
            "read1|UMI=ACGT|DUPCOUNT=1",
            "read1|UMI=ACGTACGTACGT|DUPCOUNT=0",
            "read1|UMI=ACGTACGTACGT|DUPCOUNT=x",
            "read1|DUPCOUNT=2",
        ]
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_cpm_sequence_id(value)

    def test_shm_identity_scales(self) -> None:
        self.assertAlmostEqual(shm_from_identity(98.5), 1.5)
        self.assertAlmostEqual(shm_from_identity(0.985), 1.5)
        self.assertEqual(shm_from_identity(101.0), 0.0)

    def test_canonical_modes_are_explicit(self) -> None:
        self.assertIsNone(junction_aa_rejection_reason("AR", "rg_reference"))
        self.assertEqual(
            junction_aa_rejection_reason("AR", "cpm_conservative"),
            "drop_junction_aa_length_outside_5_40",
        )
        self.assertIsNone(junction_aa_rejection_reason("CARGW", "cpm_conservative"))
        self.assertEqual(
            junction_aa_rejection_reason("AARGW", "cpm_conservative"),
            "drop_junction_aa_not_cysteine_started",
        )

    def test_family_and_supporting_read_counts_are_separate(self) -> None:
        root = Path(__file__).resolve().parent / "_runtime"
        root.mkdir(exist_ok=True)
        try:
            tsv = root / "sample.umi_exact.igblast.airr.tsv"
            rows = [
                row("r1|UMI=AAAAAAAAAAAA|DUPCOUNT=3"),
                row("r2|UMI=CCCCCCCCCCCC|DUPCOUNT=2"),
                row("r3|UMI=GGGGGGGGGGGG|DUPCOUNT=7", productive="F"),
            ]
            with tsv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            config = AnalysisConfig(
                input_path=tsv,
                output_dir=root,
                sample="sample",
                cache_path=root / "pgen_cache.tsv",
            )
            result = read_and_aggregate(config, lambda _message: None)
            (
                aa_family_counts,
                aa_supporting_counts,
                _j_to_shm,
                _j_to_aa,
                j_to_family_count,
                j_to_supporting_count,
                _j_to_vlen,
                accepted_rows,
                stats,
            ) = result

            self.assertEqual(stats["input_exact_umi_families"], 3)
            self.assertEqual(stats["input_supporting_reads"], 12)
            self.assertEqual(stats["kept_exact_umi_families"], 2)
            self.assertEqual(stats["kept_supporting_reads"], 5)
            self.assertEqual(aa_family_counts["CAR"], 2)
            self.assertEqual(aa_supporting_counts["CAR"], 5)
            self.assertEqual(sum(j_to_family_count.values()), 2)
            self.assertEqual(sum(j_to_supporting_count.values()), 5)
            self.assertEqual([item["supporting_read_count"] for item in accepted_rows], [3, 2])
            self.assertTrue(all(item["umi_family_count"] == 1 for item in accepted_rows))
        finally:
            for child in root.iterdir():
                if child.is_file():
                    child.unlink()
            root.rmdir()


if __name__ == "__main__":
    unittest.main()
