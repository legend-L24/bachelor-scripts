#!/usr/bin/python3

import json
import re
import sys
from fractions import Fraction

def bound_chk(bound):
    bound = \
        [[float(Fraction(x)) for x in s.split(",")] for s in bound.split(";")]
    bound = bound, [b if b[0] == b[1] else [0.0, 1.0] for b in bound]
    return bound[0] == bound[1]

def dic_key(k):
    ret = re.match("([0-9]+)(.*)", k[0]).groups()
    return int(ret[0]), ret[1]

def main():
    for k, v in sorted(json.load(sys.stdin).items(), key = dic_key):
        ret, v = [], v["pos"]
        for p in v[:-1]:
            if not bound_chk(p[2]):
                ret.append(p[0])
        if ret:
            sys.stdout.write("%s: %s\n" % (k, " ".join(ret)))

if __name__ == "__main__":
    main()

