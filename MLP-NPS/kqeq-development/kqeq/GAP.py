import numpy as np

from ase.io import *
from ase.data import covalent_radii
from ase.units import Bohr, Hartree
import random
from quippy.descriptors import Descriptor
from kqeq.data import atomic_numbers
from kqeq.funct import get_energies

class GAP():
    """Class for Kernel Charge Equilibration Models

    Parameters
    ----------
        Kernel: obj
            Kernel object with the function kernel_matrix
        scale_atsize: float
            Scaling factor for covalent radii, used to define atom sizes
        radius_type: string
            'rcov' or 'qeq' (default)
        sparse: bool
            Setting training to use full or sparse training set
        sparse_count: int
            Number of training points used for sparsed training
        sparse_method: string
            Method used to specify representative set "FPS", or "CUR" (default) 
    """

    def __init__(self, Kernel=None, scale_atsize=1.0, sparse=False, sparse_count = 1000, sparse_method = "CUR"):
        
        if Kernel == None:
            print("Specify kernel!")
            exit()
        # inputs and defaults
        self.scale_atsize = scale_atsize
        self.Kernel       = Kernel
        self.weights      = None
        self.sparse       = sparse
        self.K_train     = self.Kernel.kernel_matrix(kerneltype='training') #Kernelmatrix(self,self.training_set)
        if self.sparse == True:
            self.sparse_method = sparse_method
            if sparse_count == None:
                if  self.K_train.shape[0] < 1000:
                    self.sparse_count =  self.K_train.shape[0]
                    print("Training set has less than 1000 points, training will be done with the full set")
                else:
                    self.sparse_count = 1000
            else:
                self.sparse_count = int(sparse_count)
            if sparse_method == "FPS":
                self.sparse_points = self._farthest_point(self.K_train,self.sparse_count,1)
                Kernel._create_representationFPS(self.sparse_points)
            elif sparse_method == "CUR":
                self.sparse_points = self.Kernel._create_representationCUR(self.sparse_count)
            else:
                print("Specify sparse method!")
            self.K_nm, K_mm = self.Kernel.kernel_matrix(kerneltype='training')
            for a in range(K_mm.shape[0]):
                K_mm[a,a] = K_mm[a,a] + 0.0000001 
            self.K_mm = K_mm
        else:
            self.K_mm = self.K_train
            self.K_nm = self.K_train



    def train(self, lambda_reg =0.1, perAtom = False, atom_energy = None, energy_keyword = "energy"):
        L = np.zeros((len(self.Kernel.training_set_descriptors),len(self.Kernel.training_set)))
        count_c = 0
        count_r = 0
        E_ref = get_energies(mols=self.Kernel.training_set,atom_energy = atom_energy, energy_keyword = energy_keyword)
        for mol in self.Kernel.training_set:
            for at in mol:
                L[count_r,count_c] = 1
                count_r += 1
            count_c += 1
        K_nm = self.K_nm
        K_mm = self.K_mm   
        Lambda = np.linalg.inv(lambda_reg*np.eye(len(E_ref)))
        up = np.linalg.multi_dot([K_nm.T,L,Lambda,E_ref])
        down = np.linalg.multi_dot([K_nm.T, L,Lambda,L.T,K_nm]) + K_mm
        weights = np.linalg.solve(down,up)
        self.weights = weights

        

    def calculate(self, mol):
        #K = self.Kernel.kernel_matrix(mol_set1 = [mol],kerneltype="predicting")
        K,dKdr  = self.Kernel.calculate_function(mol_set1 = mol)
        energy = np.matmul(K,self.weights)
        nAt  = len(mol)
        f_k  = np.zeros((nAt,3))
        for j in range(nAt):
            for direction in range(3):
                eneg_drj = np.matmul(dKdr[j][direction],self.weights)
                for i in range(nAt):
                    f_k[j,direction] += -eneg_drj[i]
        results = {'energy':np.sum(energy)*Hartree,'forces':f_k * Hartree}
        return results

    def write_sparse_xml(self):
        with open("GAP.xml.314","w") as xml:
            for atom in self.Kernel.representing_set_descriptors:
                for el in atom:
                    xml.write(f"{str(el)}\n")
            
    def write_xml_weights(self):
        id = 1
        with open("GAP.xml","w") as xml:
            for atom in self.weights:
                xml.write(f'<sparseX i="{id}" alpha="{atom}" sparseCutoff="1.0000000000000000"/>\n')
                id += 1


    def save_model(self,name="kqeq", gap_version=1631710424, atom_energy=[]):
        sparse_len = len(self.Kernel.representing_set_descriptors)
        sparse_dim = len(self.Kernel.representing_set_descriptors[0])
        with open(f"{name}.xml.sparseX.kqeq","w") as xml:
            for atom in self.Kernel.representing_set_descriptors:
                for el in atom:
                    xml.write(f"{str(el)}\n")
        with open(f"{name}.xml","w") as xml:
            xml.write(f"<{name}>\n")
            xml.write(f'<Potential label="{name}" init_args="IP GAP label={name}"/>\n')
            xml.write(f'<GAP_params label="{name}" gap_version="{gap_version}">\n')
            xml.write(f'  <GAP_data do_core="F">\n')
            for element in atomic_numbers:
                if element not in atom_energy:
                    xml.write(f'    <e0 Z="{atomic_numbers[element]}" value=".00000000000000000E+000"/>\n')
                else:
                    xml.write(f'    <e0 Z="{atomic_numbers[element]}" value="{atom_energy[element]}"/>\n')
            xml.write(f'  </GAP_data>\n')
            xml.write(f'  <gpSparse label="{name}" n_coordinate="1">\n')
            xml.write(f'    <gpCoordinates label="{name}1" dimensions="{sparse_dim}" signal_variance="1.0" signal_mean=".00000000000000000E+000" sparsified="T" n_permutations="1" covariance_type="2" zeta="2.0000000000000000" n_sparseX="{sparse_len}" sparseX_filename="{name}.xml.sparseX.kqeq" sparseX_md5sum="2519cb796dd238f6aaf967cdf0dbb96e">\n')
            xml.write(f'            <descriptor>{self.Kernel.Descriptor_description}</descriptor>\n')
            id = 1
            for atom in self.weights:
                xml.write(f'      <sparseX i="{id}" alpha="{atom}" sparseCutoff="1.0000000000000000"/>\n')
                id += 1
            xml.write(f'    </gpCoordinates>\n')
            xml.write(f'  </gpSparse>\n')
            xml.write(f'</GAP_params>\n')
            xml.write(f'</{name}>\n')
