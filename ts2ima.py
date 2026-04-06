#!/usr/bin/env python3
"""
ts2ima.py — time x wavelength (TS2IMA analogue)

Thin wrapper around ts2dynspec.py:
  --x wavelength --y time
"""
from __future__ import annotations
import sys
from ts2dynspec import main

if __name__ == "__main__":
    sys.exit(main(["--x", "wavelength", "--y", "time", *sys.argv[1:]]))
