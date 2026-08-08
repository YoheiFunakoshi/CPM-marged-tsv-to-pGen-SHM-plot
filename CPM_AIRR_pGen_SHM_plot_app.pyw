# -*- coding: utf-8 -*-
from __future__ import annotations

import importlib
import queue
import threading
import traceback
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from cpm_airr_pgen_shm_plot import AnalysisConfig, run_analysis, safe_sample_name


APP_DIR = Path(__file__).resolve().parent


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CPM exact-UMI AIRR TSV -> pGen / SHM / pGen-SHM plot")
        self.geometry("980x720")
        self.minsize(820, 560)

        self.log_queue: queue.Queue[str] = queue.Queue()
        self.worker: threading.Thread | None = None

        self.input_var = tk.StringVar()
        self.outdir_var = tk.StringVar()
        self.sample_var = tk.StringVar()
        self.cache_var = tk.StringVar()
        self.include_supporting_var = tk.BooleanVar(value=False)
        self.recalc_pgen_var = tk.BooleanVar(value=True)
        self.pgen_workers_var = tk.StringVar(value="6")
        self.vlen_var = tk.StringVar(value="0")
        self.xlim_var = tk.StringVar(value="-30,-5")
        self.ylim_var = tk.StringVar(value="0,15")
        self.bw_var = tk.StringVar(value="0.8")
        self.canonical_var = tk.StringVar(value="RG reference: letters only, no stop/X")

        self._build_ui()
        self.after(100, self._drain_log_queue)

    def _build_ui(self):
        root = tk.Frame(self, padx=14, pady=12)
        root.pack(fill=tk.BOTH, expand=True)

        row = 0
        unit_note = tk.Label(
            root,
            text=(
                "Primary unit: exact UMI family (1 AIRR row).  "
                "DUPCOUNT is a separate supporting merged-read count."
            ),
            fg="#0b4f8a",
            anchor="w",
            justify="left",
        )
        unit_note.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))

        row += 1
        tk.Label(root, text="AIRR TSV (IgBLAST outfmt 19)").grid(row=row, column=0, sticky="w")
        tk.Entry(root, textvariable=self.input_var).grid(row=row, column=1, sticky="ew", padx=6)
        tk.Button(root, text="Browse", command=self.browse_input).grid(row=row, column=2, sticky="ew")

        row += 1
        tk.Label(root, text="Output folder").grid(row=row, column=0, sticky="w", pady=(8, 0))
        tk.Entry(root, textvariable=self.outdir_var).grid(row=row, column=1, sticky="ew", padx=6, pady=(8, 0))
        tk.Button(root, text="Browse", command=self.browse_outdir).grid(row=row, column=2, sticky="ew", pady=(8, 0))

        row += 1
        tk.Label(root, text="Sample name").grid(row=row, column=0, sticky="w", pady=(8, 0))
        tk.Entry(root, textvariable=self.sample_var).grid(row=row, column=1, sticky="ew", padx=6, pady=(8, 0))
        tk.Label(root, text="used in titles/output names").grid(row=row, column=2, sticky="w", pady=(8, 0))

        row += 1
        tk.Label(root, text="pGen cache TSV").grid(row=row, column=0, sticky="w", pady=(8, 0))
        tk.Entry(root, textvariable=self.cache_var).grid(row=row, column=1, sticky="ew", padx=6, pady=(8, 0))
        tk.Button(root, text="Default", command=self.set_default_cache).grid(row=row, column=2, sticky="ew", pady=(8, 0))

        opts = tk.LabelFrame(root, text="Options", padx=10, pady=8)
        opts.grid(row=row + 1, column=0, columnspan=3, sticky="ew", pady=(12, 0))
        opts.columnconfigure(1, weight=1)
        opts.columnconfigure(3, weight=1)
        tk.Checkbutton(
            opts,
            text="Also create DUPCOUNT supporting-read-weighted plots (secondary unit)",
            variable=self.include_supporting_var,
        ).grid(
            row=0, column=0, columnspan=4, sticky="w"
        )
        tk.Checkbutton(opts, text="Recalculate all pGen (ignore existing cache)", variable=self.recalc_pgen_var).grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(6, 0)
        )
        tk.Label(opts, text="pGen workers").grid(row=1, column=2, sticky="e", padx=(18, 4), pady=(6, 0))
        tk.Entry(opts, textvariable=self.pgen_workers_var, width=14).grid(row=1, column=3, sticky="w", pady=(6, 0))
        tk.Label(opts, text="Min V alignment length").grid(row=2, column=0, sticky="w", pady=(6, 0))
        tk.Entry(opts, textvariable=self.vlen_var, width=12).grid(row=2, column=1, sticky="w", pady=(6, 0))
        tk.Label(opts, text="KDE xlim").grid(row=2, column=2, sticky="e", padx=(18, 4), pady=(6, 0))
        tk.Entry(opts, textvariable=self.xlim_var, width=14).grid(row=2, column=3, sticky="w", pady=(6, 0))
        tk.Label(opts, text="KDE ylim").grid(row=3, column=0, sticky="w", pady=(6, 0))
        tk.Entry(opts, textvariable=self.ylim_var, width=12).grid(row=3, column=1, sticky="w", pady=(6, 0))
        tk.Label(opts, text="Bandwidth").grid(row=3, column=2, sticky="e", padx=(18, 4), pady=(6, 0))
        tk.Entry(opts, textvariable=self.bw_var, width=14).grid(row=3, column=3, sticky="w", pady=(6, 0))
        tk.Label(opts, text="Junction AA filter").grid(row=4, column=0, sticky="w", pady=(6, 0))
        ttk.Combobox(
            opts,
            textvariable=self.canonical_var,
            state="readonly",
            values=(
                "RG reference: letters only, no stop/X",
                "CPM conservative: C...F/W, length 5-40",
            ),
            width=45,
        ).grid(row=4, column=1, columnspan=3, sticky="w", pady=(6, 0))

        actions = tk.Frame(root)
        actions.grid(row=row + 2, column=0, columnspan=3, sticky="ew", pady=(12, 0))
        tk.Button(actions, text="Check setup", command=self.check_setup).pack(side=tk.LEFT)
        self.run_button = tk.Button(actions, text="Run pGen + SHM + pGen-SHM plot", command=self.run_clicked)
        self.run_button.pack(side=tk.LEFT, padx=(8, 0))

        row += 3
        tk.Label(root, text="Log").grid(row=row, column=0, sticky="w", pady=(12, 0))
        row += 1
        self.log_text = tk.Text(root, height=22, wrap="word")
        self.log_text.grid(row=row, column=0, columnspan=3, sticky="nsew", pady=(4, 0))
        scroll = tk.Scrollbar(root, command=self.log_text.yview)
        scroll.grid(row=row, column=3, sticky="ns", pady=(4, 0))
        self.log_text.configure(yscrollcommand=scroll.set)

        root.columnconfigure(1, weight=1)
        root.rowconfigure(row, weight=1)

    def browse_input(self):
        path = filedialog.askopenfilename(
            title="Select AIRR TSV",
            filetypes=[("AIRR TSV", "*.tsv *.txt *.zip"), ("All files", "*.*")],
        )
        if not path:
            return
        p = Path(path)
        self.input_var.set(str(p))
        if not self.outdir_var.get().strip():
            self.outdir_var.set(str(p.parent))
        if not self.sample_var.get().strip():
            self.sample_var.set(safe_sample_name(p))
        if not self.cache_var.get().strip():
            self.cache_var.set(str(Path(self.outdir_var.get()) / "pgen_cache.tsv"))

    def browse_outdir(self):
        path = filedialog.askdirectory(title="Select output folder")
        if path:
            self.outdir_var.set(path)
            if not self.cache_var.get().strip() or Path(self.cache_var.get()).name == "pgen_cache.tsv":
                self.cache_var.set(str(Path(path) / "pgen_cache.tsv"))

    def set_default_cache(self):
        outdir = self.outdir_var.get().strip()
        if not outdir:
            messagebox.showwarning("Output folder", "Select output folder first.")
            return
        self.cache_var.set(str(Path(outdir) / "pgen_cache.tsv"))

    def check_setup(self):
        self.clear_log()
        packages = ["numpy", "matplotlib", "scipy", "olga", "openpyxl"]
        ok = True
        for package in packages:
            try:
                module = importlib.import_module(package)
                version = getattr(module, "__version__", "OK")
                self.log(f"{package}: OK {version}")
            except Exception as exc:
                ok = False
                self.log(f"{package}: MISSING ({exc})")
        try:
            import olga

            model_dir = Path(olga.__file__).resolve().parent / "default_models" / "human_B_heavy"
            if model_dir.exists():
                self.log(f"OLGA human_B_heavy model: OK {model_dir}")
            else:
                ok = False
                self.log(f"OLGA human_B_heavy model: MISSING {model_dir}")
        except Exception:
            ok = False
        if ok:
            messagebox.showinfo("Check setup", "Setup OK.")
        else:
            messagebox.showerror("Check setup", "Some dependencies are missing. See log.")

    def validate_config(self) -> AnalysisConfig:
        input_path = Path(self.input_var.get().strip())
        if not input_path.exists():
            raise ValueError("AIRR TSV file not found.")
        outdir_text = self.outdir_var.get().strip()
        output_dir = Path(outdir_text) if outdir_text else input_path.parent
        sample = self.sample_var.get().strip() or safe_sample_name(input_path)
        cache_text = self.cache_var.get().strip()
        cache_path = Path(cache_text) if cache_text else output_dir / "pgen_cache.tsv"
        try:
            min_vlen = int(self.vlen_var.get().strip() or "0")
        except ValueError as exc:
            raise ValueError("Min V alignment length must be an integer.") from exc
        try:
            pgen_workers = int(self.pgen_workers_var.get().strip() or "6")
        except ValueError as exc:
            raise ValueError("pGen workers must be an integer.") from exc
        if pgen_workers < 1:
            raise ValueError("pGen workers must be >= 1.")
        xlim = parse_pair(self.xlim_var.get().strip(), "KDE xlim")
        ylim = parse_pair(self.ylim_var.get().strip(), "KDE ylim")
        try:
            bw = float(self.bw_var.get().strip() or "0.8")
        except ValueError as exc:
            raise ValueError("Bandwidth must be numeric.") from exc
        if bw <= 0:
            raise ValueError("Bandwidth must be > 0.")
        return AnalysisConfig(
            input_path=input_path,
            output_dir=output_dir,
            sample=sample,
            cache_path=cache_path,
            include_supporting_read_outputs=self.include_supporting_var.get(),
            recalculate_pgen=self.recalc_pgen_var.get(),
            pgen_workers=pgen_workers,
            min_v_align_len=max(0, min_vlen),
            locus="IGH",
            canonical_mode=(
                "cpm_conservative"
                if self.canonical_var.get().startswith("CPM conservative")
                else "rg_reference"
            ),
            xlim=xlim,
            ylim=ylim,
            bw_factor=bw,
            prefix=safe_sample_name(Path(sample)),
        )

    def run_clicked(self):
        if self.worker and self.worker.is_alive():
            messagebox.showwarning("Running", "Analysis is already running.")
            return
        try:
            config = self.validate_config()
        except Exception as exc:
            messagebox.showerror("Input", str(exc))
            return
        self.clear_log()
        self.run_button.configure(state=tk.DISABLED)
        self.log("Starting analysis...")

        def target():
            try:
                outputs = run_analysis(config, log=self.thread_log)
                self.thread_log("")
                self.thread_log("Complete.")
                for key, value in outputs.items():
                    self.thread_log(f"{key}: {value}")
                self.log_queue.put("__DONE__")
            except Exception:
                self.thread_log(traceback.format_exc())
                self.log_queue.put("__ERROR__")

        self.worker = threading.Thread(target=target, daemon=True)
        self.worker.start()

    def thread_log(self, message: str):
        self.log_queue.put(message)

    def _drain_log_queue(self):
        try:
            while True:
                message = self.log_queue.get_nowait()
                if message == "__DONE__":
                    self.run_button.configure(state=tk.NORMAL)
                    messagebox.showinfo("Complete", "Analysis complete.")
                elif message == "__ERROR__":
                    self.run_button.configure(state=tk.NORMAL)
                    messagebox.showerror("Error", "Analysis failed. See log.")
                else:
                    self.log(message)
        except queue.Empty:
            pass
        self.after(100, self._drain_log_queue)

    def clear_log(self):
        self.log_text.delete("1.0", tk.END)

    def log(self, message: str):
        self.log_text.insert(tk.END, str(message) + "\n")
        self.log_text.see(tk.END)
        self.update_idletasks()


def parse_pair(text: str, label: str) -> tuple[float, float]:
    try:
        left, right = [float(x.strip()) for x in text.split(",", 1)]
    except Exception as exc:
        raise ValueError(f"{label} must be two comma-separated numbers, e.g. -30,-5") from exc
    if right <= left:
        raise ValueError(f"{label}: max must be greater than min.")
    return left, right


if __name__ == "__main__":
    App().mainloop()
