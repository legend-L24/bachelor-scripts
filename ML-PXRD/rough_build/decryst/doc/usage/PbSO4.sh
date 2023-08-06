#!/bin/sh -e
#
# This example assumes that you have successfully installed decryst (including
# the `.json' files, or otherwise you will additionally need to set
# `$DECR_EXTRA'), and have the scripts from `doc/scripts' in your $PATH.  You
# should also have created a directory, entered this directory, copied
# `doc/examples/PbSO4.txt' to `pso.txt', and prepared a `hosts.conf' file in the
# same format as `doc/examples/hosts.conf'.  Additionally, if you enable remote
# computation, `decr_stat', `decr_lsa' and `decr_sas' should also be installed
# on all remote machines.
#
# Please do note that this script is not intended to be run as a whole, except
# when you want to measure the performance of decryst.  It is intended as a
# "recipe" to show you the typical workflow of decryst, which you would need to
# understand anyway.  This is not difficult: if there is a step you do not
# understand, stop and observe what this step does, and perhaps consult the
# documentation of an utility or the source code of a helper scipt.  After you
# finally understand everything, you are encouraged to read
# `doc/examples/0009563.sh'.

# Comment out the light atoms in `pso.txt'.
sed -i '/^O2-/ s/^/#/; /:2-/ s/^/#/' pso.txt
# Create a directory for solving the heavy atoms.
mkdir ps
# Dump the currently possible `cryst' files into the `ps' directory.
# Using `doc/scripts/do_dump.sh', the following can be simplified into:
# $ echo pso.txt | do_dump.sh ps
decr_utils comb pso.txt | awk '{ print "ps/pso," $0 }' |
	decr_utils dump pso.txt > ps/index.list

# Get the statistics of these `cryst' files.
# Using `doc/scripts/do_stat.sh', the following can be simplified into:
# $ do_stat.sh ps
PYTHONUNBUFFERED=1 decr_utils stat stat.sh hosts.conf < ps/index.list |
	tee ps/index.tmp | do_rank.sh 0 > ps/index.rank0
# Solve the heavy atoms now; before doing this, you can comment out lines in
# `ps/index.rank0' with too low scores, eg. with
# $ awk '{ print ($2 < 5 ? "#" : "") $0 }' < ps/index.rank0 > ps/index.tmp
# $ cp ps/index.tmp ps/index.rank0
# Using `doc/scripts/do_optim.sh', the following can be simplified into:
# $ do_optim.sh ps
grep -v '^#' < ps/index.rank0 | cut -f 2-3 |
	PYTHONUNBUFFERED=1 decr_utils optim optim.sh hosts.conf |
	tee ps/index.tmp | do_rank.sh 1 > ps/index.rank1

# Merge the results with the commented `decr' file; before doing this, you can
# comment out lines in `ps/index.rank1' with too low scores in the same way as
# with `ps/index.rank0'.
# Using `doc/scripts/do_merge.sh', the following can be simplified into:
# $ do_merge.sh pso.txt ps
do_filter.sh < ps/index.rank1 | decr_utils merge pso.txt > /dev/null
# Uncomment the light atoms in `pso.txt' and the merged `decr' files.
sed -i '/^#O2-/ s/^#//; /:2-/ s/^#//' pso.txt
do_filter.sh < ps/index.rank1 | sed 's/\.cr$/.txt/' |
	xargs -n 100 sed -i '/^#O2-/ s/^#//; /:2-/ s/^#//'
# Create a directory for solving the light atoms.
mkdir pso
# Dump the currently possible `cryst' files into the `pso' directory.
# Using `doc/scripts/do_dump.sh', the following can be simplified into:
# $ do_filter.sh < ps/index.rank1 | sed 's/\.cr$/.txt/' | do_dump.sh pso
do_filter.sh < ps/index.rank1 | sed 's/\.cr$/.txt/' | while read f; do
	decr_utils comb "$f" | awk '{
		print "pso/'"$(basename $f | sed -r 's@\.[^/.]+$@,@')"'" $0
	}' | decr_utils dump "$f" >> pso/index.list
done

# Get the statistics of these `cryst' files.
# $ do_stat.sh pso
PYTHONUNBUFFERED=1 decr_utils stat stat.sh hosts.conf < pso/index.list |
	tee pso/index.tmp | do_rank.sh 0 > pso/index.rank0
# Solve the light atoms now; before doing this, you can comment out lines in
# `pso/index.rank0' with too low scores in the same way as with
# `ps/index.rank0'.
# $ do_optim.sh pso
grep -v '^#' < pso/index.rank0 | cut -f 2-3 |
	PYTHONUNBUFFERED=1 decr_utils optim optim.sh hosts.conf |
	tee pso/index.tmp | do_rank.sh 1 > pso/index.rank1

# Merge the results with the uncommented `decr' file; before doing this, you can
# comment out lines in `pso/index.rank1' with too low scores in the same way as
# with `ps/index.rank0'.
# $ do_merge.sh pso.txt pso
do_filter.sh < pso/index.rank1 | decr_utils merge pso.txt > /dev/null
# Clean up; now you can inspect the `.cif' files in `pso', according to the
# ranks in `pso/index.rank1'.
rm ps/index.tmp pso/index.tmp

