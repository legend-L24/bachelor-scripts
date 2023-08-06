#!/bin/sh -
# Wrapper that runs `decr_mcs' on the specified machine.
# (cf. `doc/usage/cmdline.txt'.)

if [ "$2" != '127.0.0.1' ]; then
	exec ssh.sh "$2" decr_mcs "$1"
elif [ -n "$4" ]; then
	exec ssh.sh "$4" decr_mcs "$1"
else
	exec decr_mcs "$1"
fi

