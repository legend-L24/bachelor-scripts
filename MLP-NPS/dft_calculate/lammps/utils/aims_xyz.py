#!/u/ytli/anaconda3/envs/lyt/bin/python
import numpy as np
from ase.io import read, write
from ase import Atoms
from pathlib import Path

# this unit transfer is necessary in CP2K rather than in FHI-aims

eV = 2.72113838565563E+01
angstrom = 5.29177208590000E-01
GPa = 160.21766208

def cp2kout_to_aseatoms(file_name):
    files = list(open(file_name))
    if len(files)==0: return 0
    if ("Have a nice day" in files[-2]):
        chemical_symbols=[]
        positions=[]
        forces=[]
        stress_tensor = []
        for idx, line in enumerate(files):
            if "The structure contains" in line:
                total_atom_number = int(files[idx].strip("\n").split()[3])
            if "Input geometry" in line:
                pos=files[idx+7:idx+total_atom_number+7]
                cell=[files[idx+i+2].strip("\n").split()[-3:] for i in range(0,3)]
                for idx, atom_pos in enumerate(pos):
                    chemical_symbols.append(atom_pos.strip("\n").split()[-4])
                    positions.append(atom_pos.strip("\n").split()[-3:])
            if "Final atomic structure" in line:
                pos=files[idx+6:idx+total_atom_number+6]
                cell=[files[idx+i+2].strip("\n").split()[-3:] for i in range(0,3)]
                cell=np.array([[float(i) for i in j] for j in cell])
                for idx, atom_pos in enumerate(pos):
                    chemical_symbols.append(atom_pos.strip("\n").split()[-1])
                    positions.append(atom_pos.strip("\n").split()[1:4])
            if "| Total energy of the DFT / Hartree-Fock s.c.f. calculation" in line:
                energy = float(line.strip("\n").split()[-2])
            if " Total atomic forces (unitary forces cleaned)" in line:
                frc=files[idx+1:idx+total_atom_number+1]
                forces = []
                for idx, atom_pos in enumerate(frc):
                    forces.append(atom_pos.strip("\n").split()[-3:])
        positions=np.array([list(map(eval, ls)) for ls in positions])
        forces=np.array([list(map(eval, ls)) for ls in forces])
        ats=Atoms(chemical_symbols, positions=positions)
        ats.set_cell(cell)
        ats.set_array('forces', forces)
        ats.info['energy'] = energy
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
for i in range(0, 11):
    #print(i)
    str_file = Path("./output__" + str(i))
    if str_file.exists():
        try:
            A = cp2kout_to_aseatoms("output__" +str(i))
            if A != 0:
                pos.append(A)
            else: print(i)
        except:
            print("one wrong file:", i)
print(len(pos))
write("gama.xyz", pos)
