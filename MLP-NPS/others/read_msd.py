#!/u/ytli/anaconda3/envs/quip/bin/python

f = open("log.lammps","r")
read_lines = f.readlines()
msd_ls = []
check_start_1 = 0
check_start = 0
for line in read_lines:
    if "Step          Temp         c_msd[4]      v_twopoint     v_fitslope" in line:
        check_start_1 = 1
    else: check_start_1 = 0
    if check_start:
        try:
            #print(line)
            a = [float(i) for i in line.split()]
            msd_ls.append(line)
        except:
            #print("meet the end")
            check_start = 0
    if check_start_1: check_start = 1
f.close()
with open("msd.txt","w") as f:
    f.writelines(msd_ls)
