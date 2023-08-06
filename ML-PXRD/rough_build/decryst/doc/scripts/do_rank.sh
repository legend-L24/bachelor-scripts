#!/bin/sh -

# Helper that ranks the Wyckoff combinations according to some heuristics.
# Mode 1: simply sort according to the FOM factor (cf. the reference below).
# Mode 0: models with zero dof first, then sort according to the FOM factor.
#
# X.-D. Deng and C. Dong. "EPCryst: a computer program for solving crystal
# structures from powder diffraction data". Journal of Applied Crystallography,
# 2011, 44(1): 230-237.

[ "$1" = '0' ] && sopt='-k1n'
awk '{
	if ($5 == "-" || $6 == "-") print 0 " " (2 / $7) "\t" $0
	else print 1 " " (1 / $7 + 10 * $6 / $5) "\t" $0
}' | sort $sopt -k2gr

