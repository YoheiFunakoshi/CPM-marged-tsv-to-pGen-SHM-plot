@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\pythonw.exe" (
  echo Python environment not found: .venv
  echo Run setup_env.bat first.
  pause
  exit /b 1
)
".venv\Scripts\pythonw.exe" "CPM_AIRR_pGen_SHM_plot_app.pyw"
