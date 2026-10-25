"""
run_all.py -- run the whole pipeline in order (steps 01-10).

Usage:  python scripts/run_all.py
        python scripts/run_all.py --from 07     (start from a given step)

Steps are discovered by globbing NN_*.py, so adding a new numbered script
makes it part of the run automatically.

Steps 01-09 build the classified table. Step 10 is read-only: it recomputes
the headline numbers under alternative assumptions and writes to
results/robustness/, leaving data_clean/ and results/*.csv untouched.

Stops at the first step that fails. Step 02 only calls the gnomAD API for
rsIDs/variants not already saved, so with complete files it finishes quickly.
"""

import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
STEPS = sorted(p for p in SCRIPTS.glob("[0-9][0-9]_*.py"))

start = "01"
if "--from" in sys.argv:
    start = sys.argv[sys.argv.index("--from") + 1].zfill(2)

for script in STEPS:
    if script.name[:2] < start:
        continue
    print(f"\n{'=' * 70}\nRUNNING {script.name}\n{'=' * 70}", flush=True)
    result = subprocess.run([sys.executable, str(script)])
    if result.returncode != 0:
        sys.exit(f"\nFAILED at {script.name} (exit code {result.returncode})")

print("\nAll steps finished.")
