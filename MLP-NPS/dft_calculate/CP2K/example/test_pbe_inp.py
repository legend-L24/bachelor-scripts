#!/usr/bin/env python
import numpy as np
from ase.io import read, write

filename = "refer_110.xyz"

rel_max = 110+20+20
cutoff_max = 1000+200+200

rel_cut_orig = 90
cutoff_orig = 800
k_pack_ls = ["1 1 1","2 2 2","3 3 3", "4 4 4", "5 5 5"]


atom = read("./"+filename)
for jdx in range(0, 5):# last 0-8
    print("this is jdx", jdx)
    k_pack = k_pack_ls[2]
    cutoff = cutoff_orig+jdx*100
    rel_cut = rel_cut_orig+jdx*10
    if jdx == 20:
        rel_cut = rel_max
        cutoff = cutoff_max    
    cell = atom.get_cell()
    pos  = atom.get_positions()
    cs   = atom.get_chemical_symbols()

    w_str = """
@set data /data/home/ch3_105/lyt/cp2k_data
@set proj_name labeling

&GLOBAL
  PROJECT_NAME ${proj_name}
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL MEDIUM
  WALLTIME 860000
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &PRINT
      &FORCES ON
      &END FORCES
  &END PRINT
  STRESS_TENSOR ANALYTICAL
  &DFT
    BASIS_SET_FILE_NAME ${data}/BASIS_MOLOPT
    POTENTIAL_FILE_NAME ${data}/GTH_POTENTIALS
    WFN_RESTART_FILE_NAME ./${proj_name}-RESTART.wfn
    CHARGE 0
    UKS F 
    &MGRID
      CUTOFF %d
      REL_CUTOFF %d
    &END MGRID
    &KPOINTS
      SCHEME MONKHORST-PACK %s
    &END KPOINTS
    &QS
      EPS_DEFAULT 1.0E-13
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 2
    &END QS
    &SCF
         SCF_GUESS ATOMIC
         EPS_SCF 1.0E-6
         MAX_SCF 300

         ADDED_MOS 2
         &DIAGONALIZATION
            ALGORITHM STANDARD
            EPS_ADAPT 0.01
         &END DIAGONALIZATION
         &SMEAR  ON
            METHOD FERMI_DIRAC
            ELECTRONIC_TEMPERATURE [K] 300
         &END SMEAR

         &MIXING
            METHOD BROYDEN_MIXING
            ALPHA 0.2
            BETA 1.5
            NBROYDEN 8
         &END MIXING

      &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &vdW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME ${data}/dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END vdW_POTENTIAL
   &END XC
    &PRINT
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
"""%(cutoff, rel_cut, k_pack) 
    w_str += "      A %12.6f %12.6f %12.6f \n" %(cell[0][0], cell[0][1], cell[0][2])
    w_str += "      B %12.6f %12.6f %12.6f \n" %(cell[1][0], cell[1][1], cell[1][2])
    w_str += "      C %12.6f %12.6f %12.6f \n" %(cell[2][0], cell[2][1], cell[2][2])
    w_str += """
    &END CELL
    &COORD
"""
    for idx, ps in enumerate(pos):
        w_str += "%2s %12.6f %12.6f %12.6f\n" %(cs[idx], ps[0], ps[1], ps[2])

    w_str +=""" 
    &END COORD
    &KIND O
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND C
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND S
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND N
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q5
    &END KIND
    &KIND P
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q5
    &END KIND
    &KIND F
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q7
    &END KIND
    &KIND Li
        BASIS_SET DZVP-MOLOPT-SR-GTH
        POTENTIAL GTH-PBE-q3
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

"""
    f = open(filename.strip(".xyz") + "_" + str(jdx) + ".inp", "w")
    f.write(w_str)
    f.close()


