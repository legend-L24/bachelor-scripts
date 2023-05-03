import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({"font.size":12})
rcParams.update({"figure.autolayout":True})


plot_channel_volumes=True
if plot_channel_volumes:
    fig, axis = plt.subplots(4, figsize=(4,14), sharey=True)
    Ts = [400,500,600,700]
    for T, ax in zip(Ts,axis):
        li3_data1 = np.loadtxt("vol_data/vol_Li3PS4_{}.txt".format(T), skiprows=1)
        li3_data2 = np.loadtxt("vol_data/vol2_Li3PS4_{}.txt".format(T), skiprows=1)
        li3_data3 = np.loadtxt("vol_data/vol3_Li3PS4_{}.txt".format(T), skiprows=1)
        li7_data1 = np.loadtxt("vol_data/vol_Li7P3S11_{}.txt".format(T), skiprows=1)
        li7_data2 = np.loadtxt("vol_data/vol2_Li7P3S11_{}.txt".format(T), skiprows=1)
        li7_data3 = np.loadtxt("vol_data/vol3_Li7P3S11_{}.txt".format(T), skiprows=1)
        li4_data1 = np.loadtxt("vol_data/vol_Li4P2S7_{}.txt".format(T), skiprows=1)
        li4_data2 = np.loadtxt("vol_data/vol2_Li4P2S7_{}.txt".format(T), skiprows=1)
        li4_data3 = np.loadtxt("vol_data/vol3_Li4P2S7_{}.txt".format(T), skiprows=1)

        ax.plot(li3_data1[:,0],  li3_data1[:,1], color="b", label="Li3PS4")
        ax.plot(li3_data2[:,0],  li3_data2[:,1], color="b")
        ax.plot(li3_data3[:,0],  li3_data3[:,1], color="b")

        ax.plot(li7_data1[:,0], li7_data1[:,1], color="r", label="Li7PS11")
        ax.plot(li7_data2[:,0], li7_data2[:,1], color="r")
        ax.plot(li7_data3[:,0], li7_data3[:,1], color="r")

        ax.plot(li4_data1[:,0],  li4_data1[:,1],  color="g", label="Li4P2S7")
        ax.plot(li4_data2[:, 0], li4_data2[:, 1], color="g")
        ax.plot(li4_data3[:, 0], li4_data3[:, 1], color="g")

        ax.set_xlabel("threshold")
        ax.set_ylabel("relative volume")
        ax.legend()
        ax.set_title("T = {} K".format(T))
    plt.savefig("volume_analysis.pdf ")
    plt.show()

