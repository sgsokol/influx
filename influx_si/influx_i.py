#! /usr/bin/env python3
from os import path
with open(path.join(path.dirname(__file__), "influx_s.py"), "r") as f:
    exec(f.read())
