#!/bin/sh -
# Wrapper that runs `decr_sas' on the specified machine.
# (cf. `doc/usage/cmdline.txt'.)

if [ "$2" != '127.0.0.1' ]; then
	exec ssh.sh "$2" decr_sas "$1" "tcp://0.0.0.0:$3"
elif [ -n "$4" ]; then
	exec ssh.sh -L "$3:127.0.0.1:$5" "$4" decr_sas "$1" "tcp://127.0.0.1:$5"
else
	exec decr_sas "$1" "tcp://127.0.0.1:$3"
fi

