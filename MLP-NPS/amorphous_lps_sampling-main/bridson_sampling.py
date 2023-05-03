from ase.io import read,write
from grids import phosphate_grid
import random
from utils.periodic_bridson import *
from utils.fps import *


def rotate_random(bb, cm=[0, 0, 0]):
    """Rotates a component by random Euler angles
    atoms: input Atoms object
    cm: center of mass"""
    com = bb.get_center_of_mass()
    bb.positions -= com
    euler_angles = [180,90,180]
    angles = [random.random()*euler_angle for euler_angle in euler_angles]
    bb.euler_rotate(phi = angles[0], theta = angles[1], psi = angles[2], center = cm)
    return bb


def get_density(atoms):
    """Get density of cell in g/cm^3"""
    masses = atoms.get_masses()
    N_A = 6.022 * 10 ** 23
    mass = sum([mass / N_A for mass in masses])
    volume = atoms.get_volume() * 10 ** (-24)
    density = mass / volume
    return density


class BuildingBlockSet:
    def __init__(self, building_blocks):
        self.bblist = []
        for building_block_kind in building_blocks.keys():
            for i in range(building_blocks[building_block_kind]):
                self.bblist.append(read("building_blocks/{}.xyz".format(building_block_kind)))


class AnionSampler:
    def __init__(self, params):
        self.density = params["initial_density"]
        self.cell = np.array(params["cell"])
        self.cation = params["cation"]
        chargelist = {"ps4": 3, "p2s7": 4, "p2s6": 4, "p2s6_": 2}
        self.building_blocks = params["building_blocks"]
        self.building_block_set = BuildingBlockSet(self.building_blocks)
        self.n_cations = sum([chargelist[building_block]*self.building_blocks[building_block] for building_block in self.building_blocks.keys()])
        #try:
        self.axis = params["iface"]["axis"]

        self.gridpoints, self.cell = phosphate_grid(self.cell, self.building_blocks, self.n_cations,
                                                        self.cation, self.density)
        #except:
        #    print("no iface")
        #    self.gridpoints, self.cell = phosphate_grid(self.cell, self.building_blocks, self.n_cations, self.cation, self.density)
        self.building_blocks[self.cation] = self.n_cations
        self.geometry=Atoms()

    def get_gridpoint(self):
        pos = random.choice(self.gridpoints)
        idx = self.gridpoints.index(pos)
        del self.gridpoints[idx]
        return pos

    def populate_phosphates(self):
        npoints = len(self.gridpoints)
        component_list = []
        atoms = Atoms()
        self.gridpoints = np.ndarray.tolist(self.gridpoints)
        for building_block in self.building_block_set.bblist:
            pos = self.get_gridpoint()
            building_block = rotate_random(building_block)
            building_block.translate(pos)
            component_list.append(list(range(len(atoms), len(atoms) + len(building_block))))
            atoms.extend(building_block)
        atoms.set_cell(self.cell)
        atoms.set_pbc((True, True, True))
        atoms.wrap()
        self.geometry=atoms
        print('Occupied gridpoints:', npoints - len(self.gridpoints), 'empty gridpoints:', len(self.gridpoints))
        return self.geometry, component_list


def populate_phosphates(gridpoints, cell):
    """
    populate the phosphorous gridpoints with ps4 tetrahedrons, which are randomly rotated
    """
    ps4 = read("building_blocks/ps4.xyz")
    atoms = Atoms()
    gridpoints = np.ndarray.tolist(gridpoints)
    for gridpoint in gridpoints:
        building_block = rotate_random(ps4)
        building_block.translate(gridpoint)
        atoms.extend(building_block)
    atoms.set_cell(cell)
    return atoms


def sample_void(cell, cutoff):
    """
    sample the gb vac structure with periodic bridson sampling,
    select phosphorous positions with periodic fps, place rotated tetrahedrons on phosphorous sites and add
    lithium to remaining sites
    """
    data_points = Bridson_sampling_periodic(dims=np.diag(cell),
                                                                 radius=cutoff, k=30,
                                                                 initial_points=[], axis=0)
    xatoms = Atoms("X"*len(data_points), positions=data_points, cell=cell)
    mod = len(data_points) % 4
    if mod == 0:
        data_points = data_points
    else:
        # do not cut the last points, but remove points randomly
        delete_mask = []
        l = list(range(len(data_points)))
        for i in range(mod):
            choice = np.random.choice(l)
            delete_mask.append(choice)
            l.remove(choice)
        data_points = np.delete(data_points, delete_mask, axis=0)
    n_sampled = len(data_points)
    p_gridpoints, li_gridpoints = conditional_fps(data_points, n_samples=int(n_sampled / 4), initial_inds=[],
                                                  axis=0, dims=np.diag(cell))
    write("example.xyz",xatoms)
    sampled = populate_phosphates(p_gridpoints, cell)
    sampled.extend(Atoms("Na" * len(li_gridpoints), cell=cell, positions=li_gridpoints))
    sampled.pbc = [1, 1, 1]
    sampled.wrap()
    return(sampled)
    # #it_compr = IterativeCompressor(sampled, axis, vac=vac, idx=j, gapdir=self.gapdir)
    # #atoms = it_compr.iter_opt()


cell = np.array([[25, 0, 0], [0, 25, 0], [0, 0, 25]])
cutoff = 3.3
sampled = []
for i in range(1):
    sampled.append(sample_void(cell, cutoff))
write("bridson_sampled.xyz", sampled)
print(get_density(sampled[0]))
#view(sampled)

# for el1 in ["P", "Na"]:
#     for el2 in ["P", "Na"]:
#         rdf, dist = get_partial_rdf(sampled[0], 6, 50, el1, el2)
#         plt.plot(dist, rdf)
#         plt.xlabel("dist")
#         plt.ylabel("rdf")
#         plt.title(f"{el1}-{el2}")
#         plt.show()
#
#
# pot = quippy.potential.Potential(param_filename="tests/amorphous_lps_sampling_copy/gap/gp_2b_soap.xml")
# cell = np.array([[20, 0, 0], [0, 20, 0], [0, 0, 20]])
# cutoffs = np.arange(3, 4, 0.1)
# densities = []
# f_av = []
# for cutoff in cutoffs:
#     print(cutoff)
#     sampled = sample_void(cell, cutoff)
#     densities.append(get_density(sampled))
#     pot.calculate(sampled, ["forces", "energy"])
#     f_av.append(np.average(np.abs(pot.results["forces"])))
# plt.plot(cutoffs, densities)
# plt.xlabel("cutoff / A")
# plt.ylabel("density / g/cm3")
# plt.figure()
# plt.plot(densities, f_av)
# plt.xlabel("density / g/cm3")
# plt.ylabel("average force / eV/A")
# plt.show()




