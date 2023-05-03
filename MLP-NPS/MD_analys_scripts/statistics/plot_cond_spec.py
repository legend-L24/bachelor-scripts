import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.stats import pearsonr
cmap = cm.get_cmap('coolwarm')
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'figure.autolayout': True})

#print("only p2s6")
#files = ["Li3PS4.txt", "Li7P3S11.txt", "Li4P2S7.txt"]
#savenames = ["li3ps4_cond_spec.pdf", "li7p3s11_cond_spec.pdf", "li4p2s7_cond_spec.pdf"]
#markers = ["x", "o", "v"]
#for j,file in enumerate(files):
#    data = np.loadtxt(file)
#    ind = 0
#    fig = plt.figure(figsize=(5.5, 5))
#    ax = fig.add_subplot(111)
#    for i,T in  enumerate(["400", "500", "600", "700"]):
#        plt.scatter(data[ind:ind+20, 3], 0.01*data[ind:ind+20, 2], label=T+" K", color=cmap(i/4), marker="X")
#        pearson = pearsonr(data[ind:ind+20, 3], 0.01*data[ind:ind+20, 2])
#        print("Pearson correlation and p value {:.2f}  {:.2f} ".format(pearson[0], pearson[1]))
#        plt.xlabel("P$_2$S$_6^{4-}$ content / %")
#        plt.ylabel("Conductivity / S/cm")
#        ind += 20
#    ax.set_aspect(1. / ax.get_data_ratio())
#    plt.legend()
#    plt.savefig("/home/huss/MA/Report/"+savenames[j])
#    plt.figure()
#
#
#print("only p2s7")
#files = ["Li3PS4.txt", "Li4P2S7.txt", "Li7P3S11.txt"]
#files2 = ["Li3PS4_comp.txt", "Li4P2S7_comp.txt", "Li7P3S11_comp.txt"]
#savenames = ["li3ps4_cond_spec2.pdf", "li7p3s11_cond_spec2.pdf", "li4p2s7_cond_spec2.pdf"]
#markers = ["x", "o", "v"]
#for j,file in enumerate(files):
#    data = np.loadtxt(file)
#    data2 = np.loadtxt(files2[j])
#    ind = 0
#    fig = plt.figure(figsize=(5.5, 5))
#    ax = fig.add_subplot(111)
#    for i,T in  enumerate(["400", "500", "600", "700"]):
#        plt.scatter(data2[ind:ind+20, 2], 0.01*data[ind:ind+20, 2], label=T+" K", color=cmap(i/4), marker="X")
#        pearson = pearsonr(data2[ind:ind+20, 2], 0.01*data[ind:ind+20, 2])
#        print("Pearson correlation and p value {:.2f}  {:.2f} ".format(pearson[0], pearson[1]))
#        plt.xlabel("P$_2$S$_7^{4-}$ content / %")
#        plt.ylabel("Conductivity / S/cm")
#        ind += 20
#    #ax.set_aspect(1. / ax.get_data_ratio())
#    plt.legend()
#    plt.savefig("/home/huss/MA/Report/"+savenames[j])
#    plt.figure()


print("sum of p2s6 and p2s7")
files = ["Li3PS4.txt", "Li4P2S7.txt", "Li7P3S11.txt"]
files2 = ["Li3PS4_comp.txt",  "Li7P3S11_comp.txt", "Li4P2S7_comp.txt"]
savenames = ["li3ps4_sig_comp.pdf", "li7p3s11_sig_comp.pdf", "li4p2s7_sig_comp.pdf"]
markers = ["x", "o", "v"]
for j,file in enumerate(files):
    data = np.loadtxt(file)
    data2 = np.loadtxt(files2[j])
    ind = 0
    fig = plt.figure(figsize=(5.5, 5))
    ax = fig.add_subplot(111)
    for i,T in  enumerate(["400", "500", "600", "700"]):
        plt.scatter((data2[ind:ind+20, 2]+data[ind:ind+20, 3])*100, 0.01*data[ind:ind+20, 2], label=T+" K", color=cmap(i/4), marker="X")
        pearson = pearsonr(data2[ind:ind+20, 2]+data[ind:ind+20, 3], 0.01*data[ind:ind+20, 2])
        print("Pearson correlation and p value {:.2f}  {:.2f} ".format(pearson[0], pearson[1]))
        plt.xlabel("P$_2$S$_n^{4-}$ content / %")
        plt.ylabel("Conductivity / S/cm")
        ind += 20
    ax.set_aspect(1. / ax.get_data_ratio())
    if j ==0:
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [3,2,1,0]
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])
    plt.savefig("/home/huss/MA/Report/"+savenames[j])
plt.show()

# """isotropicity"""
# stds = np.zeros(20)
# msd_all = np.zeros(3)
# for i in range(20):
#     data = np.loadtxt("../plots/Li3PS4_"+str(i)+"_500_glass_all_all.txt")
#     data = data.reshape((4, int(len(data)/4))).T
#     #plt.plot(data[:,0], data[:,1])
#     #plt.plot(data[:,0], data[:,2])
#     #plt.plot(data[:,0], data[:,3])
#     #plt.figure()
#     stds[i] = np.std(data[-1, 1:])/np.mean(data[-1, 1:])
#     msd_all+=data[-1, 1:]
# print(np.mean(stds))
# print(np.std(msd_all)/np.mean(msd_all))