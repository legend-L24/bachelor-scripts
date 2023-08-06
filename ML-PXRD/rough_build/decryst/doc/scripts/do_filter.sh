#!/bin/sh -
# Helper that filters Wyckoff combinations which are commented out.
grep -v '^#' | cut -f 3

