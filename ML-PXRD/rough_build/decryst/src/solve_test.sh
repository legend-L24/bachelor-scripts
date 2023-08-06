#!/bin/sh
python3 ../python/decr_utils dump "$1" < combs.list > crysts.list
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts.conf < crysts.list > crysts.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts.conf < crysts.stat > crysts_1.optim
sort -k7g < crysts_1.optim > crysts.tmp && mv crysts.tmp crysts_1.optim
