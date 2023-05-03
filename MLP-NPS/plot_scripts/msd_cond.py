#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
from lmp_output_handler import lammps_dump_file_2_ase
from md_tools import atomic_diffusion_tools
from ase.io import read, write
from ase.data import chemical_symbols
#from utils import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt


def _load_pickle_file(filename):
    pickle_file = open(filename,'rb')
    data = pickle.load(pickle_file)
    pickle_file.close()
    return(data)


def _write_pickle_file(filename,data):
    output = open(filename, 'wb')
    pickle.dump(data, output)
    output.close()


def _fit_rates(a,b):
    """simply function doing a linear regression
    """
    fitfunc = lambda p, x: p[0] + p[1]*x # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p0 = [b[0], ((b[-1]-b[0])/(a.size))]
    pars, success = leastsq(errfunc, p0[:],args=(a,b))
    return(pars[1])


def get_diffusion_stats(time, msd):
    """ function to calculate the diffusion coefficients via Einstein relation and the
        anomalous diffusion exponent (alpha) - latter defines normal/anormal diffusion
        input:  time   = list/1-d array with time steps
                msd    = np.array of length time.size with colomns ([x^2,y^2,z^2,r^2])
        output: dfactor = diffusion coefficient
                alpha_exp = anomalous diffusion exponent
                NOTE: units are as input units (normally A^2/ps)
    """
    alpha_exp = _fit_rates(np.log((time[1:]-time[0])),np.log(msd[1:,3]))
    dfactor = (msd[-1,3]-msd[0,3])/(2*6*(time[-1]-time[0]))
    return(dfactor, alpha_exp)


def _D2conductivity_nernst(Dm,q,N,volume,temp):
    """ simple fct to calculate the conductivity
        note input is in A^2/ps output in S/m
    """
    A2cm, A2m, ps2s = 1e-08, 1e-10, 1e-12
    kB = 1.38064852e-23 #m2 kg s-2 K-1
    e = 1.60217662e-19 #el-charge
    prefactor = (e*e)/(3*volume*A2m*A2m*A2m*kB*temp)*((A2m*A2m)/ps2s)
    D_tot = sum([Dm[el]*q[el]*q[el]*N[el] for el in Dm]) #charge goes in squared!
    return(D_tot*prefactor)


def _qmsd2conductivity(time,qmsd,volume,temp):
    """ function to calculate the conductivity from the charge displacement
        after Haskins
        input:  time   = list/1-d array with time steps
                qmsd   = np.array of length time.shape
                volume = (float) in A^3
                tmep   = (float) in K
        output: sig    = conductivity in S/m
    """
#    if (time.size != -qmsavefig('{}.png'.format(filename)).size):
#        raise Exception("time and qmsd do not match")
    A2cm, A2m, ps2s = 1e-08, 1e-10, 1e-12
    kB = 1.38064852e-23 #m2 kg s-2 K-1
    e = 1.60217662e-19 #el-charge
    prefactor = (e*e)/(6*volume*A2m*A2m*A2m*kB*temp)*((A2m*A2m)/ps2s)
    sig = (qmsd[-1]-qmsd[1]) / (time[-1]-time[1]) * prefactor #output is S/m
    return(sig)


def _plot_msd_dict(filename,msd_dict,t,target_folder='output'):
    f_lin = lambda x, a, b: a + b*x
    colors=["steelblue","darkred","orange"]
    el = np.sort(list(msd_dict.keys()))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    mmax = 0
    labels = ["Li", "P", "S"]
    for i in [0]:
        msd = msd_dict[el[i]]
        mmax = max(mmax,msd[:,3].max())
        ax1.plot(t[::4],msd[::4,3],color=colors[i], label=r'%s'%chemical_symbols[el[i]])
        ax1.plot([t[0],t[-1]],[msd[1,3],msd[1,3]],color=colors[i],ls=':',lw=1) #indicating the intial offset-vibration
        np.savetxt(filename[:-3] + "_li.txt", np.array((t[::4], msd[::4, 3])))
    plt.savefig('{}.pdf'.format("li_"+filename), bbox_inches='tight')
    fig =plt.figure()
    ax1 = fig.add_subplot(111)
    for i in [1,2]:
        msd = msd_dict[el[i]]
        mmax = max(mmax,msd[:,3].max())
        ax1.plot(t[::4],msd[::4,3],color=colors[i], label=r'%s'%chemical_symbols[el[i]])
        ax1.plot([t[0],t[-1]],[msd[1,3],msd[1,3]],color=colors[i],ls=':',lw=1) #indicating the intial offset-vibration
        if i ==1:
            np.savetxt(filename[:-3]+"_p.txt", np.array((t[::4],msd[::4,3])))
        else:
            np.savetxt(filename[:-3]+"_s.txt", np.array((t[::4],msd[::4,3])))
#    ax1.yaxis.set_major_locator(MultipleLocator(0.5))
#    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax1.set_xlabel(r'lag time / ps')
    ax1.set_ylabel(r'msd / $\AA^2$')
    plt.legend(prop={"size":14}, loc = "upper left")
    plt.savefig('{}.pdf'.format("ps_"+filename), bbox_inches='tight')


def _plot_msd_all_directions(filename,msd_dict,t):
    f_lin = lambda x, a, b: a + b*x
    colors=["steelblue","darkred","orange"]
    el = ['x','y','z']#np.sort(msd_dict.keys())
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    mmax = 0

    Li_msd = msd_dict.get(3)
    x = Li_msd[:,0]
    y = Li_msd[:,1]
    z = Li_msd[:,2]
    plt.plot(t,x,label='x')
    plt.plot(t,y,label='y')
    plt.plot(t,z,label='z')
    plt.legend()
#    ax1.yaxis.set_major_locator(MultipleLocator(0.5))
#    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax1.set_xlabel(r'lag time / ps')
    ax1.set_ylabel(r'msd / $\AA^2$')
    plt.savefig('{}.pdf'.format(filename),bbox_inches='tight')
    data = np.concatenate((t,x,y,z))
    np.savetxt(filename+'_all.txt'.format(filename,), data)


def moving_av(Li_msd_dict,windowsize):
     # Also works for the(strictly invalid) cases when N is even.
    l = Li_msd_dict
    N = windowsize
    if (N//2)*2 == N:
        N = N - 1
    front = np.zeros(N//2)
    back = np.zeros(N//2)
    for i in range(1, (N//2)*2, 2):
        front[i//2] = np.convolve(l[:i], np.ones((i,))/i, mode = 'valid')
    for i in range(1, (N//2)*2, 2):
        back[i//2] = np.convolve(l[-i:], np.ones((i,))/i, mode = 'valid')
    return np.concatenate([front, np.convolve(l, np.ones((N,))/N, mode = 'valid'), back[::-1]])


def acf(x, length=200):
    return np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1]  \
        for i in range(1, length)])


def autocorr1(x,y,lags):
    '''numpy.corrcoef, partial'''

    corr=[1. if l==0 else np.corrcoef(x[l:],y[:-l])[0][1] for l in lags]
    return np.array(corr)


def disparrays(tensorpos,array2append,dispplots):
    for disp in dispplots:
        x,y,z = [],[],[]
        for i in disp:
            x.append(i[0])
            y.append(i[1])
            z.append(i[2])
        xyz = [x,y,z]
        autocor = []
        autocor = autocorr1(xyz[tensorpos[0]],xyz[tensorpos[1]],range(len(disp)-1))
        array2append.append(autocor)
    correlationmean = []
    for i in range(len(array2append[0])):
        lst = [item[i] for item in array2append]
        correlationmean.append(np.mean(lst))
    return correlationmean


def get_msd_conductivity(atoms, dumpfile, temp):
    tool = atomic_diffusion_tools() #obj-best to be inherited
    el_list = [3, 15, 16]
    # (1) msds: diffusion coefficients or conductivity (Nernst-Einstein)
    # a) raw msd
    msd_dict = tool._calc_msd_time_averaged(atoms,atypes=el_list,lag=1.0)
    msd_dict_per_atom = tool._calc_msd_per_atom(atoms,3)
    # b) time-averaged msd - smoothes out behavior (here 5% overlap shifted)
    msd_dict_av = tool._calc_msd_time_averaged(atoms,atypes=el_list,lag=1.0)
    disp_per_atoms = tool._calc_disp_per_atom(atoms,3)
    # calculate diffusion coefficients:
    # NOTE time axis: in this case I know that a snapshot is dumped every ps you
    # should adjust that automatically in your script (I have a trivial fct screening log-file for all keywords)
    # print(type(msd_dict.keys()))
    time = np.arange(msd_dict[list(msd_dict.keys())[0]][:,0].size)
    time_av = np.arange(msd_dict_av[list(msd_dict_av.keys())[0]][:,0].size)
    print("#"*20)
    dfactors, dfactors_av = {}, {} #collect diff-coeffs
    for el in el_list:
        # compute diff coeff (dfactor) via Einstein relation
        dfactor, a_exp = get_diffusion_stats(time, msd_dict[el])
        dfactor_av, a_exp_av = get_diffusion_stats(time_av, msd_dict_av[el])
        # print results per element
        print(el,":")
        print("diffusion coefficients --> raw: {:.3e}, av: {:.3e} (A^2/ps)".format(dfactor,dfactor_av))
        print("anom. diff. exponents  --> raw: {:.3e}, av: {:.3e}".format(a_exp,a_exp_av))
        dfactors.update({el:dfactor}); dfactors_av.update({el:dfactor_av})

    # calculate conductivity via Nernst-Einstein
    Nel = {el:np.where(atoms[0].get_atomic_numbers() == el)[0].size for el in el_list} # get stoichiometry
    sigma = _D2conductivity_nernst(dfactors,q={3:1,15:-2,16:-1},N=Nel,volume=atoms[0].get_volume(),temp=temp)
    sigma_av = _D2conductivity_nernst(dfactors_av,q={3:1,15:-2,16:-1},N=Nel,volume=atoms[0].get_volume(),temp=temp)
    print("\nconductivity Nernst-Einstein --> raw: {:.3e}, av: {:.3e} (S/m)\n".format(sigma,sigma_av))
    print("#"*20)
    plotting = False
    if plotting:
        _plot_msd_dict("plots/"+dumpfile+"_raw", msd_dict, time)
        _plot_msd_dict("plots/"+dumpfile+"_av", msd_dict_av, time_av)
        _plot_msd_all_directions("plots/"+dumpfile+"_all", msd_dict, time)

    calc_haskins = False
    if calc_haskins:
        ins = np.where(atoms[0].get_atomic_numbers() == 3)[0]
        msdplots = []
        dispplots = []

        for i in range(len(ins)):
            msdplots.append([])
            dispplots.append([])

        for i,n in enumerate(msd_dict_per_atom):
            for ii in range(len(ins)):
                msdplots[ii].append(n[ii][3])
#                dispplots[ii].append(disp_per_atoms[ii][3])
        for i,n in enumerate(disp_per_atoms):
           for ii in range(len(ins)):
               dispplots[ii].append(disp_per_atoms[i][ii])

        # a) raw q-msd
        qd_dict, qmsd_dict = tool._calc_current_einstein_time_averaged(atoms,charge_dict={3:1.,15:-2.,16:-1.},lag=1.0)
        # b) time-averaged msd - smoothes out behavior (here 5% overlap shifted)
        qd_dict_av, qmsd_dict_av = tool._calc_current_einstein_time_averaged(atoms,charge_dict={3:1.,15:-2.,16:-1.},lag=0.95)
        """Conductivity Haskins"""


        print("#"*20)
        for key in qmsd_dict:
            if len(qmsd_dict[key].shape) > 1:
                qmsd = qmsd_dict[key][:,3]; qmsd_av = qmsd_dict_av[key][:,3]
            else:
                qmsd = qmsd_dict[key]; qmsd_av = qmsd_dict_av[key]
            # compute conductivity (sigma) via Einstein formulation of current correlation
            sigma = _qmsd2conductivity(time,qmsd,volume=atoms[0].get_volume(),temp=temp)
            sigma_av = _qmsd2conductivity(time_av,qmsd_av,volume=atoms[0].get_volume(),temp=temp)
            # print results per element / element combination
            print(key,":")
            print("conductivity Haskins --> raw: {:.3e}, av: {:.3e} (A^2/ps)".format(sigma,sigma_av))
        print("#"*20)
    print(dfactors)
    return sigma, dfactors[3]


def conductivity_measuring_phase(dumpfile, temp):
    print("Analyse file: ", dumpfile)
    # i.e. get atoms obj from lammps dump file (here "equi_short.dump")
    # alternatively (or better) use lammpscs
    d_species = {2:15, 1:3, 3:16}; l_remove_type=[4]
    atoms = lammps_dump_file_2_ase(dumpfile,d_species,l_remove_type)
    print("Chemical composition: ", atoms[-1].get_chemical_formula())
    # for diffusion take last 5 ps, corresponds to last 20 steps
    atoms = atoms[-2500:]
#    atoms = read('cond_atoms.xyz',':')
    # alternatively load atoms file
    #atoms = read("data/example_lmp_traj_Li3OCl_1.11_T400_01/atoms_short.traj",":")

    # tool to calculate msd(s) / qmsd(s) from ase-atoms file
    tool = atomic_diffusion_tools() #obj-best to be inherited

    el_list = [3] # NOTE: this is given as an argument here i.e. which elements should be analyzed (I'd take all)

    #################################################################
    # (1) msds: diffusion coefficients or conductivity (Nernst-Einstein)
    #################################################################
    # a) raw msd
    msd_dict = tool._calc_msd_time_averaged(atoms,atypes=el_list,lag=1.0)
    msd_dict_per_atom = tool._calc_msd_per_atom(atoms,3)
    #print msd_dict
    # b) time-averaged msd - smoothes out behavior (here 5% overlap shifted)
    msd_dict_av = tool._calc_msd_time_averaged(atoms,atypes=el_list,lag=0.8)
    disp_per_atoms = tool._calc_disp_per_atom(atoms,3)
    print(disp_per_atoms.shape,disp_per_atoms[0].shape)
    # calculate diffusion coefficients:
    # NOTE time axis: in this case I know that a snapshot is dumped every ps you
    # should adjust that automatically in your script (I have a trivial fct screening log-file for all keywords)
    # print(type(msd_dict.keys()))
    time = np.arange(msd_dict[list(msd_dict.keys())[0]][:,0].size)
    time_av = np.arange(msd_dict_av[list(msd_dict_av.keys())[0]][:,0].size)
    # get diffusion coefficients for every element (el)
    print("#"*20)
    dfactors, dfactors_av = {}, {} #collect diff-coeffs
    for el in el_list:
        # compute diff coeff (dfactor) via Einstein relation
        dfactor, a_exp = get_diffusion_stats(time, msd_dict[el])
        dfactor_av, a_exp_av = get_diffusion_stats(time_av, msd_dict_av[el])
        # print results per element
        print(el,":")
        print("diffusion coefficients --> raw: {:.3e}, av: {:.3e} (A^2/ps)".format(dfactor,dfactor_av))
        print("anom. diff. exponents  --> raw: {:.3e}, av: {:.3e}".format(a_exp,a_exp_av))
        dfactors.update({el:dfactor}); dfactors_av.update({el:dfactor_av})

    # calculate conductivity via Nernst-Einstein
    Nel = {el:np.where(atoms[0].get_atomic_numbers() == el)[0].size for el in el_list} # get stoichiometry
    sigma = _D2conductivity_nernst(dfactors,q={3:1,15:-2,16:-1},N=Nel,volume=atoms[0].get_volume(),temp=temp)
    sigma_av = _D2conductivity_nernst(dfactors_av,q={3:1,15:-2,16:-1},N=Nel,volume=atoms[0].get_volume(),temp=temp)
    print("\nconductivity Nernst-Einstein --> raw: {:.3e}, av: {:.3e} (S/m)\n".format(sigma,sigma_av))
    print("#"*20)
    sig_einstein = sigma

    ins = np.where(atoms[0].get_atomic_numbers() == 3)[0]
    msdplots = []
    dispplots = []

    for i in range(len(ins)):
        msdplots.append([])
        dispplots.append([])

    for i,n in enumerate(msd_dict_per_atom):
        for ii in range(len(ins)):
            msdplots[ii].append(n[ii][3])
#            dispplots[ii].append(disp_per_atoms[ii][3])
    for i,n in enumerate(disp_per_atoms):
       for ii in range(len(ins)):
           dispplots[ii].append(disp_per_atoms[i][ii])

    ########################################################################
    # (2) calculate conductivity following the formulation by Haskins et al.
    ########################################################################
    # a) raw q-msd
    qd_dict, qmsd_dict = tool._calc_current_einstein_time_averaged(atoms,charge_dict={3:1.,15:-2.,16:-1.},lag=1.0)
    # b) time-averaged msd - smoothes out behavior (here 5% overlap shifted)
    qd_dict_av, qmsd_dict_av = tool._calc_current_einstein_time_averaged(atoms,charge_dict={3:1.,15:-2.,16:-1.},lag=0.95)

    # get conductivity after Haskins - NOTE: keywords in qmsd are to be adjusted in mdtools - in principle any ion combination
    # I refrained from that since its quite a costly computation if we talk about TB of MD data
    print("#"*20)
    sigma_haskins = []
    for key in qmsd_dict:
        # preshape depending on geometry of charge diffusion (eff. vs std.)
        if len(qmsd_dict[key].shape) > 1:
            qmsd = qmsd_dict[key][:,3]; qmsd_av = qmsd_dict_av[key][:,3]
        else:
            qmsd = qmsd_dict[key]; qmsd_av = qmsd_dict_av[key]
        # compute conductivity (sigma) via Einstein formulation of current correlation
        sigma = _qmsd2conductivity(time,qmsd,volume=atoms[0].get_volume(),temp=temp)
        sigma_haskins.append(sigma)
        sigma_av = _qmsd2conductivity(time_av,qmsd_av,volume=atoms[0].get_volume(),temp=temp)
        print(key,":")
        print("conductivity Haskins --> raw: {:.3e}, av: {:.3e} (A^2/ps)".format(sigma,sigma_av))
    print("#"*20)

    sigma_haskins = sigma_haskins[1:4]
    return sig_einstein, sigma_haskins


def read_dump(dumpfile, dtime = 200, end= -1, start = 600):
    step = int(1000/dtime)
    print("Analyse file: ",dumpfile)
    d_species = {2:15,1:3,3:16}; l_remove_type=[4]
    atoms = lammps_dump_file_2_ase(dumpfile,d_species,l_remove_type)[start:end:step]
    print("Chemical composition: ", atoms[-1].get_chemical_formula())
    return atoms


def save_conductivity(contents, file="results.txt"):
    with open(file, "a") as file:
        file.write(f' {" ".join(str(x) for x in contents)}'+"\n")


read_ensemble = True
if read_ensemble:
    for compound in ["Li3PS4", "Li4P2S7", "Li7P3S11"]:
        savefile = "conduct_data/" + compound + ".txt"
        with open(savefile, "w") as file:
            file.write("T num sigma dfactor \n")
        for T in [400, 500, 600, 700]:
            for num in range(20):
                rootdir = r"E:\MA_backup\tabea_ma_data\md_runs\production_runs\statistics_MD"
                dir = "{}/{}/{}/{}/geom.dump".format(rootdir,compound,T,num)
                print(dir)
                trj = read_dump(dir)
                sigma, dfactor = get_msd_conductivity(trj, compound+"_"+str(num)+"_"+str(T)+"_"+"glass", temp=T)
                save_conductivity([T, num, sigma, dfactor], savefile)

