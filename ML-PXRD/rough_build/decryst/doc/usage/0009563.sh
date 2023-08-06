#!/bin/sh -e
#
# If you find this example difficult to understand, try the following:
# * Simplify `doc/usage/PbSO4.sh' (as noted in the comments therein) with the
#   example helpers in `doc/scripts', keeping in mind that all these examples
#   revolve around the basic tasks summarised in `doc/usage/guide.txt'.
# * Compare the simplified PbSO4 script with this example, and (assuming you now
#   understand the lines similar to those in the simplified script) focus on
#   trying to understand what the "additional" lines do.

mkdir rgt-0
sed '/^O2-/ s/^/#/; s/^0\.25 /0.99 /' < 0009563.txt > 9563.txt
echo 9563.txt | do_dump.sh rgt-0
mkdir rgt-1
sed -i 's/^0\.99 /0.25 /' 9563.txt
echo 9563.txt | do_dump.sh rgt-1

do_stat.sh rgt-0
do_optim.sh rgt-0
mv -f rgt-0/index.rank1 rgt-0/index.tmp
awk '{ print ($9 > 0.12 ? "#" : "") $0 }' < rgt-0/index.tmp > rgt-0/index.rank1
do_filter.sh < rgt-0/index.rank1 | sed 's@.*/@rgt-1/@' > rgt-1/index.tmp
grep -Ff rgt-1/index.tmp < rgt-1/index.list > rgt-1/index.tmpp
mv -f rgt-1/index.tmpp rgt-1/index.list
do_stat.sh rgt-1
do_optim.sh rgt-1
do_merge.sh 9563.txt rgt-1
awk '{ print ($9 > 0.2 || $11 > 0.02 ? "#" : "") $0 }' \
	< rgt-1/index.rank1 > rgt-1/index.rank2

mkdir rgto-0
sed -i '/^#O2-/ s/^#//; s/^0\.25 /0.99 /' 9563.txt
do_filter.sh < rgt-1/index.rank2 | sed 's/\.cr$/.txt/' |
	xargs -n 100 sed -i '/^#O2-/ s/^#//; s/^0\.25 /0.99 /'
do_filter.sh < rgt-1/index.rank2 | sed 's/\.cr$/.txt/' | do_dump.sh rgto-0
mkdir rgto-1
sed -i 's/^0\.99 /0.25 /' 9563.txt
do_filter.sh < rgt-1/index.rank2 | sed 's/\.cr$/.txt/' |
	xargs -n 100 sed -i 's/^0\.99 /0.25 /'
do_filter.sh < rgt-1/index.rank2 | sed 's/\.cr$/.txt/' | do_dump.sh rgto-1

do_stat.sh rgto-0
mv -f rgto-0/index.rank0 rgto-0/index.tmp
awk '{ print ($9 > 0.99 ? "#" : "") $0 }' < rgto-0/index.tmp > rgto-0/index.rank0
do_optim.sh rgto-0
mv -f rgto-0/index.rank1 rgto-0/index.tmp
awk '{ print ($9 > 0.05 ? "#" : "") $0 }' < rgto-0/index.tmp > rgto-0/index.rank1
do_filter.sh < rgto-0/index.rank1 | sed 's@.*/@rgto-1/@' > rgto-1/index.tmp
grep -Ff rgto-1/index.tmp < rgto-1/index.list > rgto-1/index.tmpp
mv -f rgto-1/index.tmpp rgto-1/index.list
do_stat.sh rgto-1
do_optim.sh rgto-1
do_merge.sh 9563.txt rgto-1
awk '{ print ($11 > 0.02 ? "#" : "") $0 }' \
	< rgto-1/index.rank1 > rgto-1/index.rank2

i=0
for n in 21 11; do
	[ "$i" -eq 0 ] && c=rgto-1/index.rank2 || c=fin-$((i - 1))/index.rank1
	d=fin-"$i"
	mkdir "$d"
	do_filter.sh < "$c" | sed 's@.*/@rgto-1/@' > "$d"/index.tmp
	grep -Ff "$d"/index.tmp < rgto-1/index.rank0 |
		sed 's/rgto-1/'"$d"'/' > "$d"/index.rank0
	cp $(cat "$d"/index.tmp) "$d"
	do_optim.sh "$d"
	sed -i "$n"',$ s/^/#/' "$d"/index.rank1
	do_merge.sh 9563.txt "$d"
	i=$((i + 1))
done

c="$d"/index.rank1
d=fin-"$i"
mkdir "$d"
do_filter.sh < "$c" | sed 's@.*/@rgto-1/@' > "$d"/index.tmp
grep -Ff "$d"/index.tmp < rgto-1/index.rank0 |
	sed 's/rgto-1/'"$d"'/' > "$d"/index.list
cat "$d"/index.list | while read l; do
	name="$(echo "$l" | sed 's@.*/@@; s/\.cr$//')"
	l="$(echo "$l" | cut -f 1-2)"
	for i in $(seq 0 9); do
		cp rgto-1/"$name".cr "$d"/"$name-$i".cr
		printf '%s\t%s\n' "$l" "$d"/"$name-$i".cr >> "$d"/index.rank0
	done
done
do_optim.sh "$d"
do_merge.sh 9563.txt "$d"
awk '{ print ($11 > 0.02 ? "#" : "") $0 }' < "$d"/index.rank1 > "$d"/index.rank2

rm -f 9563.txt */index.tmp

