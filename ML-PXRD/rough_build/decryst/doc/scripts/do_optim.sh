#!/bin/sh -
# Helper that performs global optimisation on Wyckoff combinations.

grep -v '^#' < "$1"/index.rank0 | cut -f 2-3 |
	PYTHONUNBUFFERED=1 decr_utils optim optim.sh hosts.conf |
	tee "$1"/index.tmp | do_rank.sh 1 > "$1"/index.rank1

