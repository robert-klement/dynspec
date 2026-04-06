#!/usr/bin/env python3
"""ts2vima.py — velocity x time (TS2VIMA analogue).

Thin wrapper around ts2dynspec.py using explicit --mode dispatch.
This keeps one implementation of the C-parity logic while preserving the
original program entry-point name.
"""
from __future__ import annotations

import sys

from ts2dynspec import main


if __name__ == "__main__":
    # Keep positional compatibility with the historical C/MIDAS frontend wrappers:
    # this script behaves as `ts2vima` and enforces velocity/time mode.
    sys.exit(main(["--x", "velocity", "--y", "time", *sys.argv[1:]]))
