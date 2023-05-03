#!/u/ytli/anaconda3/envs/quip/bin/python
from os import mkdir
from os import popen as pp


for i in range(3,11):
    str = """#!/bin/bash
#SBATCH -o ./tjob.out.%j        #Standard output
#SBATCH -e ./tjob.err.%j        #Error output
#SBATCH -D ./                   #Initial working directory
#SBATCH -J aimsreg              #Job name
#SBATCH --nodes=1               #Nr. of compute nodes
#SBATCH --ntasks=72             #Nr. of MPI processes for the job
#SBATCH --ntasks-per-core=1     #Nr. of threads per MPI process; set always to 1
#SBATCH --mem=500000              #Memory in MB
#SBATCH --mail-type=none
#SBATCH --mail-user=vondrak@fhi-berlin.mpg.de
#SBATCH --time=24:00:00         #Wall clock time; max = 24:00:00
"""
    str += """
module purge
module load gcc/10 impi/2021.2 mkl/2021.2 fftw-mpi/3.3.9 gsl/2.4
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MKLROOT/lib/intel64"
export QUIP_ARCH=linux_x86_64_gfortran
export QUIP_ROOT=/u/ytli/apps/QUIP/quip_v0

/u/ytli/apps/QUIP/quip_v0/build/linux_x86_64_gfortran/gap_fit at_file=train.xyz sparse_jitter=1e-10 gap={distance_2b cutoff=5.0 covariance_type=ard_se cutoff_transition_width=0.6 delta=0.6 theta_uniform=1.0 sparse_method=uniform add_species=T n_sparse=20 : soap cutoff=%d cutoff_transition_width=2.0 n_max=4 l_max=12 delta=0.05 atom_sigma=0.8 zeta=2 n_sparse=200 normalise=T sparse_method=cur_points add_species covariance_type=dot_product} e0_method=isolated default_sigma={0.0005 0.02 0.1 1.0} energy_parameter_name=energy force_parameter_name=forces do_copy_at_file=F sparse_separate_file=T gp_file=gp_2b_soap.xml >& out.log
    """%(i)
    try:
        mkdir('{0}'.format(i))
    except:
        print("this file exists:", i)
    with open('{0}/gap.sbs'.format(i),"w") as f:
        f.write(str)
    pp('cp train.xyz {0} && chmod 777 {0}/gap.sbs'.format(i))
