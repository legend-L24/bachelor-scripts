#!/bin/sh -e
# See the thesis for detailed explanation.

sed -ri '/^(O2-|Na1\+)/ s/^/#/; /:2,/ s/^/#/' 158.txt
mkdir cs
echo 158.txt | do_dump.sh cs
do_stat.sh cs
do_optim.sh cs
sed -i '2,$ s/^/#/' cs/index.rank1
do_merge.sh 158.txt cs

sed -ri '/^#(O2-|Na1\+)/ s/^#//; /:2,/ s/^#//' 158.txt
do_filter.sh < cs/index.rank1 | sed 's/\.cr$/.txt/' |
	xargs -n 100 sed -ri '/^#(O2-|Na1\+)/ s/^#//; /:2,/ s/^#//'
mkdir no
do_filter.sh < cs/index.rank1 | sed 's/\.cr$/.txt/' | do_dump.sh no
do_stat.sh no
do_optim.sh no
sed -i '31,$ s/^/#/' no/index.rank1
do_merge.sh 158.txt no

mkdir csno
do_filter.sh < no/index.rank1 | sed 's@^[^,]\+,@@; s/\.cr//' |
	awk '{ print "csno/158," $0 ".cr\t" $0 }' |
	decr_utils dump 158.txt > csno/index.list
do_stat.sh csno
do_optim.sh csno
do_merge.sh 158.txt csno

mkdir fin
l="$(sed 1q < csno/index.rank1 | sed 's/csno/fin/; s/\.cr$//')"
name="$(echo "$l" | cut -f 3 | sed 's@.*/@@')"
for i in $(seq 0 9); do
	cp csno/"$name".cr fin/"$name-$i".cr
	echo "$l-$i".cr >> fin/index.rank0
done
do_optim.sh fin
do_merge.sh 158.txt fin
rm */index.tmp

