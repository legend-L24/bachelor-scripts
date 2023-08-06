#!/bin/sh
python3 ../python/decr_utils comb "$1" > combs_1.list
python3 ../python/decr_utils dump "$1" < combs_1.list > crysts_1.list
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts.conf < crysts_1.list > crysts_1.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts.conf < crysts_1.stat > crysts_1.optim
sort -k7g < crysts_1.optim > crysts_1.tmp && mv crysts_1.tmp crysts_1.optim
python3 ../python/decr_utils comb "$2" > combs_2.list
python3 ../python/decr_utils dump "$2" < combs_2.list > crysts_2.list
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts.conf < crysts_2.list > crysts_2.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts.conf < crysts_2.stat > crysts_2.optim
sort -k7g < crysts_2.optim > crysts_2.tmp && mv crysts_2.tmp crysts_2.optim
