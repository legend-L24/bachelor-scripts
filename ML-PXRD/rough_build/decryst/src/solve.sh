#!/bin/sh
python3 ../python/decr_utils comb "$1" > combs.list
python3 ../python/decr_utils dump "$1" < combs.list > crysts.list
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts_1.conf < crysts.list > crysts.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts_1.conf < crysts.stat > crysts.optim
sort -k7g < crysts.optim > crysts.tmp && mv crysts.tmp crysts.optim
python3 ../python/decr_utils merge "$1" < combs.list

