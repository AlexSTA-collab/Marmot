"""Run EdelweissFE as the analyses of the paper were run.

Two process-local settings, EdelweissFE itself untouched:
  - the correction tolerance of the scalar variables (the roller multipliers) is set to the first argument, 1.0 in
    every analysis of the paper (formerly efe_tol.py);
  - the output manager 'npydump' (double-precision nodal fields as .npy, ijes/npydump.py) is registered; only the soft
    family of the regime sweep uses it (formerly npydump/efe_npy.py).
usage: python -m ijes.efe_run <tolerance> job.inp"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from edelweissfe.config import phenomena, registry  # noqa: E402

phenomena.fieldCorrectionTolerance["scalar variables"] = float(sys.argv.pop(1))
registry._BUILTINS[("outputmanager", "npydump")] = "npydump:OutputManager"
from edelweissfe._cli._edelweissfe import main  # noqa: E402

sys.exit(main())
