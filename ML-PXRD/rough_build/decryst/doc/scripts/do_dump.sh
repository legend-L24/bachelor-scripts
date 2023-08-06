#!/bin/sh -e
# Helper that dumps Wyckoff combinations into a directory.

while read f; do
	decr_utils comb "$f" | awk '{
		print "'"$1/$(basename $f | sed -r 's@\.[^/.]+$@,@')"'" $0
	}' | decr_utils dump "$f" >> "$1"/index.list
done

