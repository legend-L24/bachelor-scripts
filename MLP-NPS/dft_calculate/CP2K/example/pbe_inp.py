#!/usr/bin/env python
import numpy as np
from ase.io import read, write

filename = "train.xyz"

atoms = read("./"+filename, "::")
for jdx, atom in enumerate(atoms):
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
    &KPOINTS
      SCHEME MONKHORST-PACK 2 2 2
    &END KPOINTS
    &MGRID
      CUTOFF 800
      REL_CUTOFF 90
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-13
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 2
    &END QS
    &SCF
      SCF_GUESS RESTART
      EPS_SCF 3.0E-7
      MAX_SCF 20
      &OUTER_SCF
        EPS_SCF 3.0E-7
        MAX_SCF 50
        #STEP_SIZE 0.03000000
      &END OUTER_SCF
      &OT
        MINIMIZER CG
        PRECONDITIONER FULL_SINGLE_INVERSE
        ENERGY_GAP 0.1
      &END OT
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
"""
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
    f = open("input_" + str(jdx) + ".inp", "w")
    f.write(w_str)
    f.close()


