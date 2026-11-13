# CoolBubble -- Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
# SPDX-License-Identifier: GPL-3.0-or-later (see the LICENSE file)
"""Re-make the plots of an existing run folder without re-running the simulation.

Usage:
    python3 python/replot.py <run folder> [<run folder> ...]
    python3 python/replot.py <run folder> --summary-only   # only the summary plots
    python3 python/replot.py <run folder> --snapshots-only # only the snapshot plots

Wall-tracking markers (r_min, r_max, gamma) are not redrawn, since they are only
known during the run.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SNAPSHOT = re.compile(r"_phi_t=(\d+)\.(dat|h5)$")


def replot(folder: Path, snapshots: bool, summary: bool) -> None:
    if not (folder / "parameters.dat").exists():
        print(f"Skipping {folder}: no parameters.dat (not a run folder?)")
        return

    if snapshots:
        files = sorted(f.name for f in folder.iterdir() if SNAPSHOT.match(f.name))
        for name in files:
            print(f"{folder.name}: {name}")
            subprocess.run([sys.executable, str(SCRIPT_DIR / "plot_FieldAndPowerSpectrum.py"),
                            "--folder", str(folder), "--filename", name], check=True)

    if summary:
        print(f"{folder.name}: summary plots")
        subprocess.run([sys.executable, str(SCRIPT_DIR / "plot_SummaryPlots.py"),
                        "--folder", str(folder)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Re-make the plots of existing run folders.")
    parser.add_argument("folders", nargs="+", type=Path, help="run folder(s) containing parameters.dat")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--summary-only", action="store_true", help="only make the summary plots")
    group.add_argument("--snapshots-only", action="store_true", help="only plot the field snapshots")
    args = parser.parse_args()

    for folder in args.folders:
        replot(folder, snapshots=not args.summary_only, summary=not args.snapshots_only)


if __name__ == "__main__":
    main()
