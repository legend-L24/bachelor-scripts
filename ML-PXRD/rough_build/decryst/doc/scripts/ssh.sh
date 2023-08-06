#!/bin/sh -

# Called by the wrappers: retry ssh(1) a few times before failing.
# When there are too many unauthenticated (yet) connections, sshd will drop new
# connections (cf. `MaxStartups' in sshd_config(5)), which is a reasonable
# security measure.  However, this means we should try several times (and use
# randomised delays to avoid too many connections retrying at the same time)
# before announcing a complete failure.

i=0
d="$(expr "$(date '+%s')" + 9)"
while true; do
	ssh "$@" && exit 0
	i="$(expr "$i" + 1)"
	[ "$(date '+%s')" -lt "$d" -a "$i" -lt 3 ] || exit 1
	sleep "$(python -c 'import random; print(1.0 + random.random())')"
done

