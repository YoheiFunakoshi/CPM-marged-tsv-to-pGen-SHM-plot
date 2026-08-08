@echo off
cd /d "%~dp0"
if "%~1"=="" (
  echo Usage: run_cli_example.bat path\to\sample.umi_exact.igblast.airr.tsv
  pause
  exit /b 1
)
".venv\Scripts\python.exe" "cpm_airr_pgen_shm_plot.py" --input "%~1" --outdir "%~dp1" --sample "%~n1"
pause
