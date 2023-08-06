#!/bin/sh -
# Helper that exports selected results in a directory.
do_filter.sh < "$2"/index.rank1 | decr_utils merge "$1" > /dev/null

