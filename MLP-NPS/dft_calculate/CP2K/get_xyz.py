#!/usr/bin/env python
import numpy as np
from ase.io import read, write
from ase import Atoms
from pathlib import Path

eV = 2.72113838565563E+01
angstrom = 5.29177208590000E-01
GPa = 160.21766208
def cp2kout_to_aseatoms(file_name):
    files = list(open(file_name))
    if "PROGRAM STOPPED IN" in files[-1]:
        chemical_symbols=[]
        positions=[]
        forces=[]
        stress_tensor = []
        for idx, line in enumerate(files):
            if "CELL| Volume" in line:
                cell=[float(files[idx+i].strip("\n").split()[-1]) for i in range(1,7)]
            if "TOTAL NUMBERS AND MAXIMUM NUMBERS" in line:
                total_atom_number = int(files[idx+3].strip("\n").split()[-1])
            if "ATOMIC COORDINATES IN angstrom" in line:
                pos=files[idx+4:idx+total_atom_number+4]
                for idx, atom_pos in enumerate(pos):
                    chemical_symbols.append(atom_pos.strip("\n").split()[2])
                    positions.append(atom_pos.strip("\n").split()[4:7])
            if "ATOMIC FORCES in [a.u.]" in line:
                frc=files[idx+3:idx+total_atom_number+3]
                for idx, atom_pos in enumerate(frc):
                    forces.append(atom_pos.strip("\n").split()[3:7])
            if "ENERGY|" in line:
                energy = float(line.strip("\n").split()[-1])
            if "STRESS TENSOR [GPa]" in line:
                strs = files[idx+3: idx+3+3]
                for idx, stress in enumerate(strs):
                    stress_tensor.append(stress.strip("\n").split()[1:4])

        positions=np.array([list(map(eval, ls)) for ls in positions])
        forces=np.array([list(map(eval, ls)) for ls in forces])
        stress_tensor=np.array([list(map(eval, ls)) for ls in stress_tensor]) # GPa
        ats=Atoms(chemical_symbols, positions=positions, cell=cell)
        ats.set_array('forces', forces * eV/angstrom)
        ats.info['energy'] = energy * eV
        ats.pbc = True
        #ats.info['stress'] = stress_tensor
        #volume = ats.get_volume()
        #virial = stress_tensor * volume / GPa
        #print(virial)
        #ats.info['virial'] = virial
        #print(stress_tensor)
        return ats
    else:
        return 0

pos =[]
for i in range(0, 300):
    #print(i)
    str_file = Path("./output_" + str(i))
    if str_file.exists():
        A = cp2kout_to_aseatoms("output_"+str(i))
        if A != 0:
            pos.append(A)
write("2.xyz", pos)
