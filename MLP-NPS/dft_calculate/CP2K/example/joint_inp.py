#!/usr/bin/env python

file_1 = "./add_data/nvt_data_0.xyz"
file_2 = "./add_data/npt_data_0.xyz"
file_name = "test.xyz"

data = []
with open(file_1,"r") as f:
    f_r = f.read()
    data.append(f_r)
with open(file_2,"r") as f:
    f_r = f.read()
    data.append(f_r)
f = open(file_name,"w")
s = ""
f.write(s.join(data))
