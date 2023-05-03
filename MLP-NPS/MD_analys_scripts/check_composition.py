from ase import atoms
from lmp_output_handler import lammps_dump_file_2_ase


def get_density(atoms):
    """Get density of cell"""
    """In the end density should 1.8g/cm^3"""
    masses = atoms.get_masses()
    N_A = 6.022*10**23
    masses = [mass/N_A for mass in masses]
    total_mass = sum(masses)
    volume = atoms.get_volume()
    volume = volume*10**(-24)
    density = total_mass/volume
    return density


d_species = {2:15,1:3,3:16}; l_remove_type=[4]
# Li10P4S15 = lammps_dump_file_2_ase('Li10P4S15/equi_1.dump',d_species,l_remove_type)
# Li10P4S15_2 = lammps_dump_file_2_ase('Li10P4S15/equi_2.dump',d_species,l_remove_type)
# print(len(Li10P4S15), len(Li10P4S15_2))
# Li11P5S18 = lammps_dump_file_2_ase('Li11P5S18/equi_Li11_1,8.dump', d_species,l_remove_type)[-1]
# Li11P5S18_2 = lammps_dump_file_2_ase('Li11P5S18/equi_Li11_1,9.dump', d_species,l_remove_type)[-1]
# # print("Li10P4S15:")
# # print(Li10P4S15.get_chemical_formula())
# # print(Li10P4S15_2.get_chemical_formula())
# # print(get_density(Li10P4S15), get_density(Li10P4S15_2))
# print("Li11P5S18:")
# print(Li11P5S18.get_chemical_formula())
# print(Li11P5S18_2.get_chemical_formula())
# print(get_density(Li11P5S18), get_density(Li11P5S18_2))

Li10P4S15_2 = lammps_dump_file_2_ase('dumps/Li10P4S15_1,9.dump', d_species,l_remove_type)[-1]
print(get_density(Li10P4S15_2))


