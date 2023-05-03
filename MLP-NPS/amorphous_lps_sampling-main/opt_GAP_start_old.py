from ase.io import read,write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
import quippy
import numpy as np
from ase import Atoms
from os import mkdir,system
from ase.constraints import UnitCellFilter
from joblib import Parallel, delayed
from ase.visualize import view
import json


class IterativeCompressor:
    def __init__(self, geom, params):
        self.geom=geom
        self.cation = params["cation"]
        self.opt_params = params["opt_params"]
        self.rho_target = params["rho_target"]
        self.gap = params["gap"]

    def get_density(self):
        """Get density of cell in g/cm^3"""
        masses = self.geom.get_masses()
        N_A = 6.022 * 10 ** 23
        mass = sum([mass / N_A for mass in masses])
        volume = self.geom.get_volume() * 10 ** (-24)
        density = mass / volume
        return density

    def get_component_list(self):
        """get a list of components to enable shrinking later without having to search for the tetrahedrons
        list is in the following form: [[0,1,2,3,4],[5]], element in atoms object 0-4 belong together
        create a string in the following form: PSSSPSSSLLP"""
        str = ""
        for el in self.geom:
            if el.symbol == self.cation:
                str = str + self.cation[0]
            else:
                str = str + el.symbol
        component_list = []
        start = 0
        for substring in ["PPSSSSSSS", "PPSSSSSS", "PSSSS", self.cation[0]]:
            component_list += [list(range(i, i + len(substring))) for i in range(start, len(str)) if
                               str.startswith(substring, i)]
            if substring == "PPSSSSSSS":
                start = component_list[-1][-1]
            if substring == "PPSSSSSS":
                start = component_list[-1][-1]
        return component_list

    def shrink_grid(self):
        """shrink cell vectors by factor, shrink also the atoms but not the phosphates
        these are moved as entire units"""
        component_list = self.get_component_list()
        cell = np.array(self.geom.cell)
        new_cell = self.opt_params["compr_factor"] * cell
        new_atoms = Atoms(cell=new_cell, pbc=(True, True, True))
        for indizes in component_list:
            Molecule = Atoms()
            for index in indizes:
                Molecule.append(self.geom[index])
            com = Molecule.get_center_of_mass()
            coeffs = np.linalg.solve(cell, com)
            new_com = np.dot(new_cell, coeffs)
            Molecule.translate(new_com - com)
            new_atoms.extend(Molecule)
        new_atoms.wrap()
        self.geom = new_atoms

    def iter_opt(self, geo_ind):
        counter = 1
        pot = quippy.potential.Potential("IP GAP label={}".format(self.gap["label_file"]),
                                         param_filename=self.gap["gapdir"])
        print("Current density: ", self.get_density())
        while self.get_density() < self.rho_target:
            if counter == 1:
                self.geom.set_calculator(pot)
                ucf = UnitCellFilter(self.geom)
                optimizer = BFGS(ucf)
                optimizer.run(fmax=self.opt_params["fmax"], steps=self.opt_params["nsteps"])
                counter += 1
                continue
            self.shrink_grid()
            print("current density: ", self.get_density())
            self.geom.set_calculator(pot)
            optimizer = BFGS(self.geom)
            optimizer.run(fmax=self.opt_params["fmax"], steps=100)
            counter += 1
        print("Run successful, Reached density: ", self.get_density())
        write("{}.xyz".format(geo_ind), self.geom)


"""read parameters from hypers file"""
geofile = "Full_Geo.xyz"
jsonpath = "params/li3ps4.json"
if __name__ == "__main__":
    Voro_file = read(geofile, ':')
    with open(jsonpath, 'r') as f:
        params = json.load(f)
    densrange = np.linspace(params["final_densrange"][0], params["final_densrange"][1], len(Voro_file))
    compressors = []
    for i, n in enumerate(Voro_file):
        params["rho_target"] = densrange[i]
        compressors.append(IterativeCompressor(geom=n, params=params))
    processed_list = Parallel(n_jobs=8)(delayed(compressor.iter_opt)(i) for i, compressor in enumerate(compressors))
    final_geom = [read("{}.xyz".format(i)) for i in range(len(compressors))]
    write("final_compressed.xyz", final_geom)
