import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
temperatures = [1000, 1200, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
files = ["plots/"+str(temperatures[i])+"_.txt" for i in range(len(temperatures))]
print(files)
for i,file in enumerate(files):
    data = np.loadtxt(file)
    plt.plot(data[0], data[1], label=str(temperatures[i])+ " K")
    plt.ylabel("msd / $\AA^2$")
    plt.xlabel("t / ps")
plt.legend()
plt.savefig("mds_vs_T.pdf")
plt.show()

temperatures = [1000, 1200, 1500]
files = ["plots/"+str(temperatures[i])+"_.txt" for i in range(len(temperatures))]
for i,file in enumerate(files):
    data = np.loadtxt(file)
    plt.plot(data[0], data[1], label=str(temperatures[i])+ " K")
    plt.ylabel("msd / $\AA^2$")
    plt.xlabel("t / ps")
plt.legend()
#plt.savefig("mds_vs_T.pdf")
plt.show()