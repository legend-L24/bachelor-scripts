import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from scipy.optimize import fmin
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.autolayout': True})
#data = np.loadtxt("results.txt")


def chi_2(params, x, y, sigy):
    m, c = params
    return sum(((y - m * x - c) / sigy) ** 2)


def err_linreg(x, y, dy, sigs = False):
    data_in=(x,y,dy)
    params0=[1,0]
    q=fmin(chi_2, params0, args=data_in)
    """error"""
    if sigs:
        N= len(x)
        D = 1/(N-1)* (sum([ (q[0]*x[i]+ q[1]- y[i])**2 for i in range(N)])) / ( N*sum([x[i]**2 for i in range(N)]) - sum([x[i]**2 for i in range(N)]))
        sig_a = np.sqrt(N*D)
        sig_b = np.sqrt(D*sum([x[i]**2 for i in range(N)]))
        return q, [sig_a, sig_b]
    else:
        return q


def extr_rt_sigma(regr, T=300):
    return 0.01*np.exp(1/T*regr[0]+regr[1])/T


def get_rt_sig_error(regr, regr_err, T= 300):
    dy_a = 0.01/(T**2)*np.exp(regr[0]/T + regr[1])
    dy_b = 0.01/(T)*np.exp(regr[0]/T + regr[1])
    return np.sqrt((dy_a*regr_err[0])**2+ (dy_b*regr_err[1])**2)


def plot_arrhenius(T, sig_array, label, ensemble = True, color="b", rt=False, tlim=False):
    if ensemble:
        sigma_av = [np.mean(i) for i in sig_array]
        std = np.array([np.std(i) for i in sigma]) * 1 / (sigma_av)
        ln_s = np.log(sigma_av * T)
        eb1 =plt.errorbar(1 / T, ln_s, yerr=std, marker="o", label=label, lw=2, capsize=3, capthick=2, linestyle="dotted", color=color)
        plt.plot(1 / T, ln_s, marker="o", color=color, lw=0)
        regr, errors = err_linreg(1 / T, ln_s, std, sigs=True)
    else:
        ln_s = np.log(sig_array*T)
        plt.plot(1/T, ln_s, "o", label=label, color=color, lw=2, linestyle="dotted")
        if len(T)>2:
            regr, errors = np.polyfit(1/T, ln_s, deg=1, cov=True)
        else:
            regr = np.polyfit(1/T, ln_s, deg=1)
        if len(T) > 2:
            errors = np.sqrt(np.diag(errors))
        else:
            errors= [0,0]
    print("eakt and rt conductivity [{:.3e},  {:.3e}]".format( -regr[0]*constants.k/constants.e, extr_rt_sigma(regr)))
    print("errors of eakt and rt conductivity [{:.3e},  {:.3e}]".format(errors[0]*constants.k/constants.e,  get_rt_sig_error(regr, errors)))
    T = np.concatenate(([300], T))
    if tlim == 1400:
        T = np.concatenate((T, [1400]))
        plt.xlim(1 / (max(T)) - .0001, 1 / (min(T) - 10))

    plt.plot(1/T, 1/T * regr[0] + regr[1], color=color)
    plt.xlabel("1/T / 1/K")
    plt.ylabel("ln($\sigma$T)")
    if rt:
        plt.plot(10 * [1 / 300], np.linspace(min(1 / T * regr[0] + regr[1]) - .5, max(ln_s) + .5, 10), color="black",
                 linestyle="dotted")
        plt.text(y =5.5, x =1/(330), s="T$_{RT}$")
    #plt.show()
'''
glass = True
if glass:
    #ax = fig.add_subplot(111)
    fig, ax = plt.subplots(figsize=(6,5))
    """li3ps4 glass"""
    data = np.loadtxt("statistics/Li3PS4.txt")
    T = np.array([400, 500, 600, 700])
    sigma = np.reshape(data[:80,2], (4, 20))
    plot_arrhenius(T, sigma, label="Li$_3$PS$_4$", color="mediumturquoise", rt=True)

    """li7p3s11"""
    data = np.loadtxt("statistics/Li7P3S11.txt")
    sigma = np.reshape(data[:80,2], (4,20))
    plot_arrhenius(T, sigma, label="Li$_7$P$_3$S$_{11}$", color="g")

    """li10p4s5"""
    data = np.loadtxt("statistics/Li4P2S7.txt")
    sigma = np.reshape(data[:80,2], (4,20))
    plot_arrhenius(T, sigma, label="Li$_{4}$P$_2$S$_{7}$", color="gold", rt=True)
    plt.savefig("arrhenius_glasses.png")
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    plt.legend()
    ax.set_aspect(1. / ax.get_data_ratio())
    plt.savefig("/home/huss/MA/Report/arrhenius_glass.pdf")

plt.figure()



"""crystals"""
print("crystals")
crystals = True
if crystals:
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    T = np.loadtxt("li7p3s11.txt")[:,0]
    sig_li7 = np.loadtxt("li7p3s11.txt")[:,2]
    plot_arrhenius(T,sig_li7, label="Li$_7$P$_3$S$_{11}$", ensemble=False, color="g", rt=True)
    plt.legend(loc="upper right")
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax.set_aspect(1. / ax.get_data_ratio())
    plt.savefig("/home/huss/MA/Report/arrhenius_li7p3s11.pdf")
    plt.figure()
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    """alpha"""
    T = np.loadtxt("alpha_true.txt")[:, 0]
    sig_alpha = np.loadtxt("alpha_true.txt")[:, 2]
    plot_arrhenius(T, sig_alpha, label="$\\alpha$-Li$_3$PS$_4$", ensemble=False, color="slateblue", tlim=1400)
    """beta"""
    T = np.loadtxt("beta_corr_lowT.txt")[:,0]
    sig_beta = np.loadtxt("beta_corr_lowT.txt")[:,2]
    plot_arrhenius(T[:3], sig_beta[:3], label="$\\beta$-Li$_3$PS$_4$", ensemble=False, color="darkblue")
    plt.plot(1/T, np.log(T*sig_beta), marker="o", color="darkblue", lw=0)
    T = np.loadtxt("gamma_cr.txt")[:,0]
    sig_gamma = np.loadtxt("gamma_cr.txt")[:,2]
    plot_arrhenius(T[:3],sig_gamma[:3], label="$\\gamma$-Li$_3$PS$_4$", ensemble=False, color="darkviolet", rt=True)
    plt.plot(1/T, np.log(T*sig_gamma), marker="o", color="darkviolet", lw=0)

    plt.legend(loc="lower left")
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax.set_aspect(1. / ax.get_data_ratio())
    plt.savefig("/home/huss/MA/Report/arrhenius_li3ps4.pdf")
plt.show()
'''
dic = "media/yli/LaCie/"


sig = [[0.6929291222312189, 0.6929291222312189], [1.0390180601894277, 1.0390180601894277], [2.3797028099694892, 2.3797028099694892], [4.5769443541284245, 4.5769443541284245], [5.232347712193188, 5.232347712193188]]
sig = [[1.3822912440387745, 1.3822912440387745], [3.6511171858454654, 3.6511171858454654], [7.476036141446458, 7.476036141446458], [10.72021170315438, 10.72021170315438], [17.876239437175855, 17.876239437175855]]
#sig = [[0.45467182433640524, 0.45467182433640524], [0.7986557674735817, 0.7986557674735817], [0.7590189389854207, 0.7590189389854207], [2.5507268842472777, 2.5507268842472777], [4.77784813349663, 4.77784813349663]]
sig = [[1.5336848477640086, 1.5336848477640086], [6.752540294015052, 6.752540294015052], [7.584222495491056, 7.584222495491056], [15.068945442082686, 15.068945442082686], [18.417076146196628, 18.417076146196628]]
all_sig = [
    [[1.5336848477640086, 1.5336848477640086], [6.752540294015052, 6.752540294015052], [7.584222495491056, 7.584222495491056], [15.068945442082686, 15.068945442082686], [18.417076146196628, 18.417076146196628]], 
    [[1.5336848477640086, 1.5336848477640086], [6.752540294015052, 6.752540294015052], [7.584222495491056, 7.584222495491056], [15.068945442082686, 15.068945442082686], [18.417076146196628, 18.417076146196628]], 
    [[0.5719718370754314, 0.5719718370754314], [1.8857350019310686, 1.8857350019310686], [5.677669825669291, 5.677669825669291], [8.290230707823259, 8.290230707823259], [17.701865154936716, 17.701865154936716]]]
all_sig = [[[1.8181346674425811, 1.8181346674425811], [4.117450783916478, 4.117450783916478], [11.900804476404225, 11.900804476404225], [14.502684433811798, 14.502684433811798], [25.755464825637628, 25.755464825637628]]]

all_sig = [
    [[1.3822912440387745, 1.3822912440387745], [3.6511171858454654, 3.6511171858454654], [7.476036141446458, 7.476036141446458], [10.72021170315438, 10.72021170315438], [17.876239437175855, 17.876239437175855]],
    [[1.5336848477640086, 1.5336848477640086], [6.752540294015052, 6.752540294015052], [7.584222495491056, 7.584222495491056], [15.068945442082686, 15.068945442082686], [18.417076146196628, 18.417076146196628]], 
    [[0.5719718370754314, 0.5719718370754314], [1.8857350019310686, 1.8857350019310686], [5.677669825669291, 5.677669825669291], [8.290230707823259, 8.290230707823259], [17.701865154936716, 17.701865154936716]]]
#all_sig = [[[1.8181346674425811, 1.8181346674425811], [4.117450783916478, 4.117450783916478], [11.900804476404225, 11.900804476404225], [14.502684433811798, 14.502684433811798], [25.755464825637628, 25.755464825637628]]]

for sig in all_sig:
    t_ls = np.array([400,500,600,700,800])
    sig_ls = np.array([i[0] for i in sig])
    plot_arrhenius(t_ls,sig_ls, label="Li$_7$P$_3$S$_{11}$", ensemble=False, color="g", rt=True)
plt.savefig("./na7p3s11.png")
'''
sig = [[0.45467182433640524, 0.45467182433640524], [0.7986557674735817, 0.7986557674735817], [0.7590189389854207, 0.7590189389854207], [2.5507268842472777, 2.5507268842472777], [4.77784813349663, 4.77784813349663]]
t_ls = np.array([400,500,600,700,800])
sig_ls = np.array([i[0] for i in sig])
plot_arrhenius(t_ls,sig_ls, label="Li$_7$P$_3$S$_{11}$", ensemble=False, color="g", rt=True)
plt.savefig("./glass_3_4.png")
'''