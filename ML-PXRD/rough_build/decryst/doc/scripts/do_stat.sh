#!/bin/sh -
# Helper that performs statistics analyses on Wyckoff combinations.

PYTHONUNBUFFERED=1 decr_utils stat stat.sh hosts.conf < "$1"/index.list |
	tee "$1"/index.tmp | do_rank.sh 0 > "$1"/index.rank0

