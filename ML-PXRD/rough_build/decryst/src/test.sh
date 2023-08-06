#!/bin/sh
python3 ../python/decr_utils stat stat.sh ../doc/examples/hosts.conf < crysts.list > crysts.stat
python3 ../python/decr_utils optim optim.sh ../doc/examples/hosts.conf < crysts.stat > crysts.optim
sort -k7g < crysts.optim > crysts.tmp && mv crysts.tmp crysts.optim
