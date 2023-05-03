from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import Calculator


class Kqeq(Calculator):
    "Calculator for Kernel Charge Equilibration Methods (KQEq)."
    def __init__(self, kqeq, **kwargs):
        """
        Parameters
        ----------
        kqeq : kernel_qeq-object
            An instance of the kernel_qeq-class for which
            training has been performed (i.e. weights are available).
        """
        super().__init__()

        # Checks
        ## we need a functional, i.e. already trained kqeq instance
        assert kqeq.weights is not None, 'A trained kqeq is required'

        # Definitions
        self.kqeq = kqeq
        self.implemented_properties = ['energy', 'forces', 'dipole_vector', 'charges']


    def calculate(self, atoms=None, properties=['energy', 'forces', 'dipole_vector', 'charges'], system_changes=all_changes):
        "Calculates properties for given Atoms-object and fills self.results."
        # basic Calculator class asks us to call its calculate-function before our implementation
        # (will set the atoms attribute)
        super().calculate(atoms, properties, system_changes)

        self.results.update(self.kqeq.calculate(atoms))

    def get_dipole_vector(self, atoms=None):
        "Return value of property ``dipole_vector`` (which is not a standard ASE property)."
        return self.get_property('dipole_vector', atoms)
