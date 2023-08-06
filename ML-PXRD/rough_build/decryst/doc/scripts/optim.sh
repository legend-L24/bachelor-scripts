#!/bin/sh -
# Wrapper that runs `decr_lsa' on the specified machine.
# (cf. `doc/usage/cmdline.txt'.)

if [ "$5" != '127.0.0.1' ]; then
	#echo "this is different 1"
        exec ssh.sh "$5" decr_lsa "$1" "$2" "$3" "$4"
elif [ -n "$7" ]; then
        #echo "this is different 2"
	exec ssh.sh "$7" decr_lsa "$1" "$2" "$3" "$4"
else
        #echo "this is different 3"
	exec decr_lsa "$1" "$2" "$3" "$4"
fi

