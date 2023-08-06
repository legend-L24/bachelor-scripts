#!/bin/sh -e
#
# Usage (following the order in the thesis):
#   /path/to/thesis.sh 0000220.txt 220
#   /path/to/thesis.sh 0000589.txt 589
#   /path/to/thesis.sh 0000428.txt 428
#   /path/to/thesis.sh 0000219.txt 219 71

mkdir $2
echo $1 | do_dump.sh $2
do_stat.sh $2
[ -z "$3" ] || sed -i "$3"',$ s/^/#/' $2/index.rank0
do_optim.sh $2
do_merge.sh $1 $2
rm -f $2/index.tmp

