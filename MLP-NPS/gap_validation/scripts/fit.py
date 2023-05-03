#!/home/ytli/anaconda3/envs/quip/bin/python
import os
os.system("/home/ytli/softwares/QUIP/build/linux_x86_64_gfortran_openmp/gap_fit default_sigma={0.006 0.033 10.0 0.0} atoms_filename=train.xyz gp_file=gp_2b_soap.xml do_copy_at_file=F \
force_parameter_name=forces energy_parameter_name=energy \
sparse_separate_file=T at_file=train.xyz verbosity=VERBOSE \
gap={distance_2b cutoff=5.5 covariance_type=ard_se add_species=T delta=.6 theta_uniform=1 n_sparse=20 sparse_method=uniform : \
\
soap_turbo l_max=3 alpha_max={{8 8 8}} atom_sigma_r={{0.5 0.5 0.5}} atom_sigma_t={{0.5 0.5 0.5}} atom_sigma_r_scaling={{0. 0. 0.}} \
atom_sigma_t_scaling={{0. 0. 0.}} zeta=4 rcut_soft=5. rcut_hard=6. basis=poly3gauss scaling_mode=polynomial amplitude_scaling={{1.0 1.0 1.0}} \
n_species=3 species_Z={{11 15 16}} radial_enhancement=1 central_index=1 central_weight={{1.0 1.0 1.0}} add_species=F \
n_sparse=2000  delta=0.1 f0=0.0 covariance_type=dot_product sparse_method=cur_points : \
\
soap_turbo l_max=3 alpha_max={{8 8 8}} atom_sigma_r={{0.5 0.5 0.5}} atom_sigma_t={{0.5 0.5 0.5}} \
atom_sigma_r_scaling={{0. 0. 0.}} atom_sigma_t_scaling={{0. 0. 0.}} zeta=4 rcut_soft=5. rcut_hard=6. \
basis=poly3gauss scaling_mode=polynomial amplitude_scaling={{1.0 1.0 1.0}} n_species=3 \
species_Z={{11 15 16}} radial_enhancement=1 central_index=2 central_weight={{1.0 1.0 1.0}} add_species=F \
n_sparse=2000 delta=0.1 f0=0.0 covariance_type=dot_product sparse_method=cur_points : \
\
soap_turbo l_max=3 alpha_max={{8 8 8}} atom_sigma_r={{0.5 0.5 0.5}} atom_sigma_t={{0.5 0.5 0.5}} \
atom_sigma_r_scaling={{0. 0. 0.}} atom_sigma_t_scaling={{0. 0. 0.}} zeta=4 rcut_soft=5. rcut_hard=6. \
basis=poly3gauss scaling_mode=polynomial amplitude_scaling={{1.0 1.0 1.0}} n_species=3 \
species_Z={{11 15 16}} radial_enhancement=1 central_index=3 central_weight={{1.0 1.0 1.0}} add_species=F \
n_sparse=2000 delta=0.1 f0=0.0 covariance_type=dot_product sparse_method=cur_points}")
