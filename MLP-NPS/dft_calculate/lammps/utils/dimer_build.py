from wfl.configset import OutputSpec
from wfl.generate import atoms_and_dimers


outputspec = OutputSpec(files='ps_dimers.xyz')

outputs = atoms_and_dimers.prepare(outputs=outputspec, atomic_numbers=[15, 16], bond_lengths={15:1.95, 16:2.01}, do_isolated_atoms=False,)
