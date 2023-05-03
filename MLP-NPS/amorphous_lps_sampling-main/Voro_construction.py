import numpy as np
from ase import Atoms, Atom
from ase.io import read,write
from ase import Atoms
from ase.visualize import view
import numpy.linalg as la
import numpy as np
import random
import matplotlib.pyplot as plt
from grids import phosphate_grid
from scipy.spatial import Voronoi,KDTree
import random
from sklearn.cluster import KMeans
import json


def rotate_random(bb, cm = [0,0,0]):
    """Rotates a component by random Euler angles
    atoms: input Atoms object
    cm: center of mass"""
    com = bb.get_center_of_mass()
    bb.positions -= com
    euler_angles = [180,90,180]
    angles = [random.random()*euler_angle for euler_angle in euler_angles]
    bb.euler_rotate(phi = angles[0], theta = angles[1], psi = angles[2], center = cm)
    return bb


class BuildingBlockSet:
    def __init__(self, building_blocks):
        self.bblist = []
        for building_block_kind in building_blocks.keys():
            for i in range(building_blocks[building_block_kind]):
                self.bblist.append(read("building_blocks/{}.xyz".format(building_block_kind)))


class Anion_Sampler:
    def __init__(self, jsonpath):
        with open(jsonpath, 'r') as f:
            params = json.load(f)
        self.density = params["initial_density"]
        self.cell = np.array(params["cell"])
        self.cation = params["cation"]
        chargelist = {"ps4": 3, "p2s7": 4, "p2s6": 4, "p2s6_": 2}
        self.building_blocks = params["building_blocks"]
        self.building_block_set = BuildingBlockSet(self.building_blocks)
        self.n_cations = sum([chargelist[building_block]*self.building_blocks[building_block] for building_block in self.building_blocks.keys()])
        self.gridpoints, self.cell = phosphate_grid(self.cell, self.building_blocks, self.n_cations, self.cation, self.density)
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


class VoronoiTesselation:
    def __init__(self, jsonpath, anion_geometry):
        with open(jsonpath, 'r') as f:
            params = json.load(f)
        self.cation = params["cation"]
        self.geometry=anion_geometry

    def cell_bounds(self):
        transposed_cell = self.geometry.cell.transpose()
        inside = []
        for i in transposed_cell:
            max_i = max(i)
            min_i = min(i)
            axes = [min_i, max_i]
            inside.append(axes)
        return inside

    def write_Anion_file(self, Anion_symbols, write_file=False):
        Anion_positions = []
        An_sym = []
        for i, n in enumerate(self.geometry.get_chemical_symbols()):
            if n in Anion_symbols:
                An_sym.append(n)
                Anion_positions.append(self.geometry.positions[i])
        Anions = Atoms('H{}'.format(len(An_sym)))
        Anions.set_chemical_symbols(An_sym)
        Anions.set_cell(self.geometry.get_cell())
        Anions.set_positions(Anion_positions)
        if write_file:
            write('Anion_grid.xyz', Anions)
        else:
            self.tesselation_grid = Anions

    def Anion_Voronoi(self, buffer_region=0):
        vor = Voronoi(self.tesselation_grid.positions)
        Voronoi_vertices = vor.vertices
        cell_boundaries = self.cell_bounds()
        inside_indices = []
        for vert in Voronoi_vertices:
            if cell_boundaries[0][0] <= vert[0] <= cell_boundaries[0][1] + buffer_region and cell_boundaries[1][
                0] <= vert[1] <= cell_boundaries[1][1] + buffer_region and cell_boundaries[2][0] <= vert[2] <= \
                    cell_boundaries[2][1] + buffer_region:
                inside_indices.append(vert)
        self.vor_vertices=inside_indices

    def Voronoi_transfer(self):
        Voronoi_atoms = Atoms('X{}'.format(len(self.vor_vertices)))
        Voronoi_atoms.set_positions(self.vor_vertices)
        self.geometry.extend(Voronoi_atoms)
        self.geometry.wrap()

    def Build_Voronoi_geometry(self, filter_dist=2.2):
        Vor_symbols = self.geometry.get_chemical_symbols()
        positions = self.geometry.get_positions()
        PPos, XPos = [], []
        for i in range(len(Vor_symbols)):
            if Vor_symbols[i] == 'P':
                PPos.append(positions[i])
            elif Vor_symbols[i] == 'X':
                XPos.append(positions[i])
        XPos = np.array(XPos)
        X_not_in = []
        for i, ii in enumerate(PPos):
            for j, n in enumerate(XPos):
                dist_XP = self.distance_3d(ii, n)
                if dist_XP < filter_dist:
                    X_not_in.append(j)
        X_not_in = np.unique(X_not_in)
        new_Xpos = []
        for i, n in enumerate(XPos):
            if i not in X_not_in:
                new_Xpos.append(n)
        del self.geometry[[atom.index for atom in self.geometry if atom.symbol == 'X']]
        XX = Atoms('X{}'.format(len(new_Xpos)))
        XX.set_positions(new_Xpos)
        self.geometry += XX

    def construct_Voronoi_cluster(self, number_of_cations):
        positions = self.geometry.get_positions()
        chemsyms = self.geometry.get_chemical_symbols()
        newpos = []
        for i, n in enumerate(chemsyms):
            if n == 'X':
                newpos.append(positions[i])
        kmeans = KMeans(n_clusters=number_of_cations)
        kmeans = kmeans.fit(newpos)
        centers = kmeans.cluster_centers_
        kdtree = KDTree(newpos)
        Lidist, Xpoints = kdtree.query(centers, 1)
        Lipoints = [newpos[ii] for ii in Xpoints]
        Li_atoms = Atoms('{}{}'.format(self.cation, number_of_cations))
        Li_atoms.set_cell(self.geometry.get_cell())
        Li_atoms.set_positions(Lipoints)
        self.geometry += Li_atoms
        del self.geometry[[atom.index for atom in self.geometry if atom.symbol == 'X']]

    def Find_nearest_Li_ion(self, Voronoi_neighbour):
        Li_neighbour = []
        for i, n in enumerate(Voronoi_neighbour):
            Li_neighbour.append(np.argmin(n))
        return Li_neighbour

    def Li_Voronoi_neighbour_list(self, Voronoi_geometry):
        Voro, Li = [], []
        Voro_neigh_Li = []
        symbs = Voronoi_geometry.get_chemical_symbols()
        for i, n in enumerate(symbs):
            if n == self.cation:
                Li.append(i)
            elif n == 'X':
                Voro.append(i)
        neighboursdists = Voronoi_geometry.get_all_distances()
        for i in Voro:
            Voro_neigh_Li.append([neighboursdists[i][j] for j in Li])
        Li_neighbour = self.Find_nearest_Li_ion(Voro_neigh_Li)
        Li_Vor_list = []
        for i in np.unique(Li_neighbour):
            Li_Vor_list.append([])
        for i, n in enumerate(Li_neighbour):
            Li_Vor_list[n].append(Voro[i])
        return Li_Vor_list

    def distance_3d(self, p1,p2):
        p1=np.array(p1)
        p2=np.array(p2)
        dist = np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2] - p2[2])**2)
        return dist


"""construction loop"""
jsonpath="params/li3ps4.json"
geoms = []
for i in range(3):
    sampler = Anion_Sampler(jsonpath)
    sampler.populate_phosphates()
    Voro=VoronoiTesselation(jsonpath, sampler.geometry)
    Voro.write_Anion_file('S')
    Voro.Anion_Voronoi()
    Voro.Voronoi_transfer()
    Voro.Build_Voronoi_geometry()
    Voro.construct_Voronoi_cluster(sampler.building_blocks[Voro.cation])
    geoms.append(Voro.geometry)
    print(Voro.geometry.symbols)
view(geoms)
write("Full_Geo.xyz", geoms)