#!/bin/sh
python3 ../python/decr_utils comb "$1" > combs_3.list
python3 ../python/decr_utils dump "$1" < combs_3.list > crysts_3.list
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts_1.conf < crysts_3.list > crysts_3.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts_1.conf < crysts_3.stat > crysts_3.optim
sort -k7g < crysts_3.optim > crysts_3.tmp && mv crysts_3.tmp crysts_3.optim
