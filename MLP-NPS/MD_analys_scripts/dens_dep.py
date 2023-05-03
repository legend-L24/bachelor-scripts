#from example_get_msd_conductivity import *
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({"figure.autolayout":True})
plt.rcParams.update({'text.usetex': True})

#evaluate_li3 = False
#if evaluate_li3:
#    T = 500
#    rho = [1.7, 1.8, 1.9]
#    for i,dens in enumerate(["rho_17", "", "rho_19"]):
#            for file in [0,1,2]:
#                dir = "../../../production_runs/glass_arrhenius/" + dens + "/li3ps4/500k/"+str(file)+"/geom.dump"
#                trj = read_dump(dir)
#                sigma, sigma_av = get_msd_conductivity(trj, "li3ps4_dens_500k_"+str(rho[i]), temp=T)
#                print(sigma)
#                save_conductivity([rho[i], file, sigma_av], file="li3ps4_dens_500k.txt")
#
#evaluate_li7= False
#if evaluate_li7:
#    T = 500
#    rho = [1.7, 1.8, 1.9]
#    for i,dens in enumerate(["rho_17", "", "rho_19"]):
#            for file in [0,1,2]:
#                dir = "../../../production_runs/glass_arrhenius/" + dens + "/li7p3s11/500k/"+str(file)+"/geom.dump"
#                trj = read_dump(dir)
#                sigma, sigma_av = get_msd_conductivity(trj, "li7p3s11_dens_500k_"+str(rho[i]), temp=T)
#                print(sigma)
#                save_conductivity([rho[i], file, sigma_av], file="density_dep/li7p3s11_dens_500k.txt")
#
plot = True
if plot:
    #matplotlib.use('TkAgg')
    fig, ax = plt.subplots(figsize=(5, 5))
    data = np.loadtxt("density_dep/li3ps4_dens_500k.txt")
    cond = data[:, 2].reshape((3,3)).T
    plt.errorbar([1.7, 1.8, 1.9], 0.01*np.mean(cond, axis=0), 0.01*np.std(cond, axis=0), label="Li$_3$PS$_4$", capsize=7, color="blue")
    data2 = np.loadtxt("density_dep/li7p3s11_dens_500k.txt")
    cond = data2[:, 2].reshape((3, 3)).T
    plt.errorbar([1.7, 1.8, 1.9], 0.01* np.mean(cond, axis=0), 0.01*np.std(cond, axis=0), label="Li$_7$P$_3$S$_{11}$", capsize=7, color="green")
    plt.xlabel("$\\rho$ / g/cm$^3$")
    plt.ylabel("$\sigma$ / S/cm")
    plt.ylim((0.0, 0.045))
    plt.legend()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_aspect(1. / ax.get_data_ratio())
    plt.savefig('dens_dep_cond.png', bbox_inches='tight', dpi=300)
    #plt.savefig("/home/huss/MA/Report/dens_dep_cond.pdf")
    plt.show()
