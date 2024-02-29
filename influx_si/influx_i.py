#! /usr/bin/env python3
from os import path
from sys import argv
d=path.dirname(argv[0])
with open(path.join(d, "influx_s.py"), "r") as f:
    exec(f.read())
