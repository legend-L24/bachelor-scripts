#!/usr/bin/python3

import json
import sys
from fractions import Fraction
from itertools import permutations
from decr_utils import pos_split, pos_expand, aff_scan, aff_fmt

def pos_fmt(pos):
    return "; ".join(",".join(aff_fmt(a) for a in p) for p in pos)

def pos_sig0(p):
    return [1 if any(a[:3]) else a[3] for a in p]

def pos_sig1(p):
    p, d = [aff_fmt(a) for a in p], {}
    for i, a in enumerate(p):
         d.setdefault(a, []).append(i)
    return sorted(v for v in d.values() if len(v) > 1)

def bound_proc(bound):
    ret = [bound.split(";"), [], Fraction(1)]
    for i, b in enumerate(ret[0]):
        b = b.split(",")
        if b[0] == b[1]:
            continue
        b = [Fraction(x) for x in b]
        ret[1].append(i)
        ret[2] *= b[1] - b[0]
    return ret

def wyck_proc(sg, w):
    pos = sg[1][w]
    pos, bret = pos_expand(sg[2], pos[0], [
        [aff_scan(a) for a in p] for p in pos_split(pos[2])
    ]), bound_proc(pos[1])
    for sf in [pos_sig0, pos_sig1]:
        sig = sf(pos[0])
        pos = [pos[0]] + [p for p in pos[1:] if sf(p) == sig]
        if sf == pos_sig0:
            pos = [[a for a in p if any(a[:3])] for p in pos]
        else:
            for i in reversed([i for v in sig for i in v[1:]]):
                for p in pos:
                    p.pop(i)
    return [pos] + bret

def pos_perm(pos, var):
    return sorted([[
        [a[i] for i in var] + [0] * (3 - len(var)) + [a[3]]
    for a in p] for p in pos])

def entry_key(e):
    k = e[0].split("; ")
    return [
        len(k), len(k[0].split(",")),
        len([x for x in "xyz" if x in k[0]]), k
    ]

def main(f):
    sgs = json.load(open(f))
    poss = {}
    for l in sys.stdin:
        sg, ws = l.split(":")
        sg, ws = [sg, sgs[sg]], ws.split()
        sg = [sg[0], dict((p[0], p[1:]) for p in sg[1]["pos"]), [
            [aff_scan(a) for a in p] for p in pos_split(sg[1]["orig"])
        ]]
        for w in ws:
            pos, bound, var, vol = wyck_proc(sg, w)
            vs = [list(v) for v in permutations(var)]
            for v in vs:
                p = pos_fmt(pos_perm(pos, v))
                if p in poss:
                    break
            else:
                v = vs[0]
                p = pos_fmt(pos_perm(pos, v))
                poss[p] = [{}, []]
            bound = [bound[i] for i in v]
            poss[p][0].setdefault(";".join(bound), vol)
            poss[p][1].append([sg[0], w, v])
    for p, [b, l] in sorted(poss.items(), key = entry_key):
        sys.stdout.write("%s\t%s\t%s\n" % ("|".join(
            e[0] for e in sorted(b.items(), key = lambda e: [-e[1], e[0]])
        ), p, "; ".join("%s,%s,%s" % (
            x[0], x[1], "".join("%u" % i for i in x[2])
        ) for x in l)))

if __name__ == "__main__":
    main(*sys.argv[1:])

