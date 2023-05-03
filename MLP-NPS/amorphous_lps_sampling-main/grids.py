import numpy as np


def phosphate_grid(cell, building_blocks, n_cations, cation, target_rho):
    masses = {"Li": 6.94, "Na": 22.99, "p2s7": 286.40, "ps4": 159.23, "p2s6":254.34, "p2s6_":254.34}
    a = cell[0]
    b = cell[1]
    c = cell[2]
    # cell volume and current density when sampled
    V = np.dot(np.cross(a, b), c)
    N_A = 6.022 * 10 ** 23
    total_mass = sum([masses[key] * building_blocks[key] / N_A for key in building_blocks.keys()])
    print(total_mass)
    total_mass += masses[cation]*n_cations/N_A
    rho = (total_mass / V) * 10 ** 24
    # rescale the cell
    cell = np.array([vect * (rho / target_rho) ** (1 / 3) for vect in [a, b, c]])
    n_points = int(2*sum(building_blocks.values()))
    a_new = cell[0]
    b_new = cell[1]
    c_new = cell[2]
    final_rho = (total_mass / np.dot(np.cross(a_new, b_new), c_new)) * 10 ** 24
    print("Sampling density is: ", final_rho, "g/cm^3")
    # rough estimate for the number of points in the cell
    # to yield a good approximate to "npoints" floor is used twice and ceil once
    n_a = int(np.floor(n_points ** (1 / 3)))
    n_b = int(np.ceil(n_points ** (1 / 3)))
    n_c = round(n_points/(n_a*n_b))
    # span 3D array with new number of gridpoints
    unit_a = a_new / n_a
    unit_b = b_new / n_b
    unit_c = c_new / n_c
    zero_point = 0.5*(unit_a+unit_b+unit_c)
    array = []
    for i in range(n_a):
        for j in range(n_b):
            for k in range(n_c):
                array.append(zero_point+i*unit_a + j*unit_b+k*unit_c)
    array = np.array(array)
    return array, cell


