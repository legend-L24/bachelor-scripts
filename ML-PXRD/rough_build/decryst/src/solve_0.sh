#!/bin/sh
python3 ../python/decr_utils comb "$1" > combs_0.list
python3 ../python/decr_utils dump "$1" < combs_0.list > crysts_0.list
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts.conf < crysts_0.list > crysts_0.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts.conf < crysts_0.stat > crysts_0.optim
sort -k7g < crysts_0.optim > crysts_0.tmp && mv crysts_0.tmp crysts_0.optim
python3 ../python/decr_utils merge "$1" < combs_0.list
