import sys
from ase.io import read, write

atoms_ls = read(sys.argv[1],"::")

def get_density(atoms):
    """Get density of cell in g/cm^3"""
    masses = atoms.geom.get_masses()
    N_A = 6.022 * 10 ** 23
    mass = sum([mass / N_A for mass in masses])
    volume = atoms.geom.get_volume() * 10 ** (-24)
    density = mass / volume
    return density
for atoms in atoms_ls:
  print(get_density(atoms)
