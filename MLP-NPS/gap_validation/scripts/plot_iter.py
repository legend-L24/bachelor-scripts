#!/home/ytli/anaconda3/envs/quip/bin/python
import matplotlib.pyplot as plt
import numpy as np

with open("error_crystal.txt", "r") as f:
	ls = f.readlines()[1:]
ls = [[float(i) for i in l.split()] for l in ls]
print(len(ls), len(ls[0]))
print(ls)
print(np.matrix(ls).shape)
ls = np.matrix(ls)
num = ls.shape[0]
train_err = ls[:,0]*1000
val_err1 = ls[:,2]*1000

lines = []
style = ['-','--', '-.']

plt.xlabel("generation")
plt.ylabel("RMSE energy/ meV/atom")

lines+=plt.plot(range(num), train_err, style[0], color='black')
lines+=plt.plot(range(num), val_err1, style[1], color='black')


#plt.axis('equal')

plt.legend(lines, ["training error", "validation error 1"], frameon=False)

plt.savefig("crystal_energy.png")

plt.clf()

train_err = ls[:,1]
val_err1 = ls[:,3]

lines = []
style = ['-','--', '-.']

plt.xlabel("generation")
plt.ylabel("RMSE force/ eV/A")

lines+=plt.plot(range(num), train_err, style[0], color='black')
lines+=plt.plot(range(num), val_err1, style[1], color='black')

#plt.axis('equal')

plt.legend(lines, ["training error", "validation error 1"], frameon=False)

plt.savefig("crystal_force.png")
