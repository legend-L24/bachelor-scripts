import numpy as np
from ase.io import read
import CG as cg 
from scipy import fftpack, ndimage
import sys
import matplotlib.pyplot as plt
#from animate_imshow import make_movie
from ase.visualize.plot import plot_atoms
from ase.visualize import view
import copy

def get_3d_density(atoms,elements,res,alpha=None,numfact=None,num=None):
    cell  = atoms.cell
    start = [0.,0.,0.]
    stop  = [cell[i,i] for i in range(len(cell))]
    res   = res

    if alpha == None:
        alpha = 0.25 #0.1 #0.25
    else:
        alpha = alpha
    if num is None:
        if numfact == None:
            num   = ((np.floor(stop)-np.round(start))/res+1).tolist()
            num   = [int(n) for n in num]
        else:
            num   = ((np.floor(stop)-np.round(start))/res+1).tolist()
            num   = [int(n)*numfact for n in num]
    
    centers = [atom.position for atom in atoms if atom.symbol in elements]

    g = cg.Grid(start=start,stop=stop,num=num)
    s = cg.Stencil(g)
    s.Gaussian(alpha)
    g = g.cyclicShape(s)
    c = np.zeros(g.shape, dtype=np.float32)
    s.place(c, list(centers))
    g,c_el = g.foldBack(c, pbcDims=None, dtype=np.float32)
    return c_el

def roll_z(atoms,density,slice_thickness):
    """slice_thickness in layers of former grid"""
    z = atoms.cell[2,2]
    #print(z)
    z_dens = np.linspace(0.,z,density.shape[2])
    z_mean = [(z_dens[i]+z_dens[i+1])/2 for i in range(len(z_dens)-1)]

    dens_xy = []
    thick = slice_thickness 
    for i in range(density.shape[2]-thick):
        dens_layer = np.sum(density[:,:,i:i+thick+1],axis=-1)
        dens_xy.append(dens_layer)
    for i in range(density.shape[2]-thick,density.shape[2]):
        dens_top = np.sum(density[:,:,i:i+thick],axis=-1)
        dens_bot = np.sum(density[:,:,:density.shape[2]-i],axis = -1)
        dens_layer = np.add(dens_top,dens_bot)
        dens_xy.append(dens_layer)

    return dens_xy


def fourier(density_xy):
    fftdata = []
    fftimage = []
    for layer in density_xy: 
        dens = np.sum(layer)
        fft2 = np.fft.fft2(layer)
        fft2 = np.fft.fftshift(fft2)
        inverse = np.fft.ifft2(fft2)
        fft2 = abs(fft2)
        fft2 = fft2/dens
        norm = np.max(fft2)
        fft2 = fft2/norm
        fftimage.append(fft2)
        fftdata.append(np.sum(fft2))
    return fftdata, fftimage


def array_2dfourier(density_xy):
    fftslab = []
    for layer in density_xy:
        dens = np.sum(layer)
        fft2 = np.fft.fft2(layer)
        fft2 = np.fft.fftshift(fft2)
        inverse = np.fft.ifft2(fft2)
        fft2 = abs(fft2)
        fft2 = fft2/dens
        norm = np.max(fft2)
        fft2 = fft2/norm
        fftslab.append(fft2)
    fftslab = np.array(fftslab)
    return fftslab

def get_binned_amorphization(atoms,elementlist, resolution,binsize,slice_thickness):
    #binned amorphization in z direction
    rho_frame = get_3d_density(atoms,elementlist,res=resolution)
    dens_xy = roll_z(atoms,rho_frame,slice_thickness)
    fft_data, fftimage = fourier(dens_xy)
    fftimage=array_2dfourier(dens_xy)

    binned = []
    for i in range(int((len(fft_data)-binsize)/binsize)):
        val = np.average(fft_data[i*binsize:(i+1)*binsize])
        if val!=val:
            val=0

        binned.append(val)
    return binned, dens_xy, fftimage


if __name__=='__main__':
    #read in desired atomic structure
    if len(sys.argv) < 2:
    	filename = raw_input("Please enter filename: ")
    else:
    	filename= sys.argv[1]
    
    atoms=read(filename)

    elementlist=['Zr', 'Y', 'La', 'Mn', 'Sr']
    amorph_profile, dens_xy, fftimage  = get_binned_amorphization(atoms, elementlist,resolution=0.1,binsize=10,slice_thickness=25)
#    print(np.shape(dens_xy))
#    print(np.shape(fftimage))
    plt.plot(range(len(amorph_profile)),amorph_profile)

#    plt.imshow(fftimage[10,:,:])
    plt.show()


