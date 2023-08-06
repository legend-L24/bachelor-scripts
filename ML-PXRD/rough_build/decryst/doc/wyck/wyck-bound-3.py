#!/usr/bin/python3

import re
import sys

def proc(sub, l):
    m = re.match(r'[^"]+"([^"]+)":.*', l)
    if m:
        sub["sg"] = m.group(1)
        return l
    m = re.match(r'([^"]+")([^"])(", [01], ")([^"]+)(".*)', l)
    if not m:
        return l
    m = list(m.groups())
    m[3] = m[3].split(";")
    e = sub.get((sub["sg"], m[1]))
    if not e:
        return l
    for i in range(3):
        if i not in e[0]:
            b = m[3][i].split(",")
            if b[0] != b[1]:
                m[3][i] = "0,1"
    b = [m[3][i] for i in e[0]]
    if b not in e[1]:
        for i, s in zip(e[0], e[1][0]):
            m[3][i] = s
    return "".join(m[:3] + [";".join(m[3])] + m[4:]) + "\n"

def main(f):
    sub = {}
    for l in sys.stdin:
        l = l.strip().split("\t")
        l = [[b.split(";") for b in l[0].split("|")], l[2].split("; ")]
        for e in l[1]:
            e = e.split(",")
            sub[(e[0], e[1])] = [[int(i) for i in e[2]], l[0]]
    for l in open(f):
        sys.stdout.write(proc(sub, l))

if __name__ == "__main__":
    main(*sys.argv[1:])

