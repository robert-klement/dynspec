#!/usr/bin/env python3
"""
ts2phima.py — phase x wavelength (TS2PHIMA analogue)

Thin wrapper around ts2dynspec.py:
  --x wavelength --y phase
"""
from __future__ import annotations
import sys
from ts2dynspec import main

if __name__ == "__main__":
    sys.exit(main(["--x", "wavelength", "--y", "phase", *sys.argv[1:]]))
