#!/usr/bin/env python3
"""
ts2vphima.py — phase x velocity (TS2VPHIMA analogue)

Thin wrapper around ts2dynspec.py:
  --x velocity --y phase
"""
from __future__ import annotations
import sys
from ts2dynspec import main

if __name__ == "__main__":
    sys.exit(main(["--x", "velocity", "--y", "phase", *sys.argv[1:]]))
