#!/usr/bin/env python3
"""
sparkto_b_raw.py — Self-contained CLI runner / alias for the 3D spherical top-hat
SPARC baryonic density & TAL factor audit (sparc_tophat.py, V19 DOI: 10.5281/zenodo.23497128).

If invoked without --data / --repo / --out arguments, defaults to the local package paths:
  --data ./data --repo . --out ./output_tophat
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import sparc_tophat

if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--data" not in argv:
        argv = ["--data", str(HERE / "data")] + argv
    if "--repo" not in argv:
        argv = ["--repo", str(HERE)] + argv
    if "--out" not in argv:
        argv = ["--out", str(HERE / "output_tophat")] + argv
    sys.argv = [sys.argv[0]] + argv
    sparc_tophat.main()
