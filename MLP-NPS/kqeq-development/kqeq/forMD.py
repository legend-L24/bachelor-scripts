import numpy as np
from ase.io import *
from ase.data import covalent_radii
from ase.units import Bohr,Hartree
from kqeq.qeq import charge_eq
from kqeq.funct import _block_diag_rect,_block_diag,_get_R, get_energies, get_charges

from sklearn.metrics import mean_squared_error

class kQEqMD():
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

    def __init__(self, Kernel=None, scale_atsize=1.0, calculate_kernel_forces=True, 
                include_LJ=False, radius_type='qeq', sparse=False, hard_lib=None):
        
        if Kernel == None:
            print("Specify kernel!")
            exit()
        # inputs and defaults
        self.scale_atsize = scale_atsize
        self.Kernel       = Kernel
        self.calculate_kernel_forces = calculate_kernel_forces
        self.LJ = include_LJ
        self.radius_type  = radius_type
        self.weights      = None
        self.hard_lib     = hard_lib
        self.sparse       = sparse
        if self.sparse == True:
            self.Kernel.sparse = True
        
    def predict(self,predict_set,predict_system_charges=None):
        """Predicts dipole vectors, charges and electronegativities for list of atoms objects.
 
        Parameters
        ----------
            predict_set : list
                List of atoms objects
            predict_system_charges : list
                List of integers representing charges of system (if None, charges are set up to be zero)

        Returns
        -------
            dipole: array
                1D array of dipole vector elements
            charge: array
                1D array of charges
            eneg: array
                1D array of electronegativities

        """ 
        if predict_system_charges == None:
            predict_system_charges = [0 for temp in predict_set]   
        K = self.Kernel.kernel_matrix(mol_set1 = predict_set,kerneltype="predicting")
        A, X, O, R, Xback = self._build_prediction(predict_set, predict_system_charges)
        eneg  = np.matmul(K,self.weights)
        eneg_tot = (np.matmul(X,np.transpose(eneg)))
        eneg_tot = eneg_tot + O
        charge_temp = np.matmul(A,-eneg_tot)
        charge = np.matmul(Xback,charge_temp)
    
        dipole  = np.matmul(R,charge_temp)
        return dipole,charge,eneg



    def calculate(self, mol, charge=0):
        """Calculates charges, energy, forces and dipole_vector for single atoms object.

        Parameters
        ----------
            mol : obj
                Atoms object
            charge : int
                Charge of mol

        Returns
        -------
            results: dict
                Dictionary of results

        """

        K,dKdr  = self.Kernel.calculate_function(mol_set1 = mol)#,predict_set=self.training_set,kerneltype="predicting") 
        eneg   = np.matmul(K,self.weights)
        qe     = charge_eq(mol,Qtot=charge,eneg=eneg,scale_atsize=self.scale_atsize,DoLJ=self.LJ,radius_type=self.radius_type, hard=self.hard_lib)
        qe.calc_Charges()
        qe.calc_Eqeq()
        charges  = qe.q
        #if self.LJ:
        #    qe.calc_ELJ()
        mu = np.matmul(_get_R(mol),qe.q)
        
        energy  = qe.E 
        forces  = qe.f *(Hartree/Bohr)
        nAt  = len(charges)
        f_k  = np.zeros((nAt,3))
        for j in range(nAt):
            for direction in range(3):
                eneg_drj = np.matmul(dKdr[j][direction],self.weights)
                for i in range(nAt):
                    f_k[j,direction] += -charges[i]*eneg_drj[i]
        forces += f_k * Hartree#/Bohr
        results = {'charges':charges,'energy':energy*Hartree,'dipole_vector':mu,'forces':forces}
        return results

    
    def calculate_old(self, mol, charge=0):
        """Calculates charges, energy, forces and dipole_vector for single atoms object.

        Parameters
        ----------
            mol : obj
                Atoms object
            charge : int
                Charge of mol

        Returns
        -------
            results: dict
                Dictionary of results

        """
        K = self.Kernel.kernel_matrix(mol_set1 = [mol],kerneltype="predicting")#,predict_set=self.training_set,kerneltype="predicting") 
        eneg   = np.matmul(K,self.weights)
        
        qe     = charge_eq(mol,Qtot=charge,eneg=eneg,scale_atsize=self.scale_atsize,DoLJ=self.LJ,radius_type=self.radius_type, hard=self.hard_lib)
        qe.calc_Charges()
        qe.calc_Eqeq()
        charges  = qe.q
        if self.LJ:
            qe.calc_ELJ()
            
        mu = np.matmul(_get_R(mol),qe.q)
        energy  = qe.E 
        forces  = qe.f *(Hartree/Bohr)
        if self.calculate_kernel_forces:
            forces += self._kernel_forces(mol,charges) * Hartree#/Bohr
        results = {'charges':charges,'energy':energy*Hartree,'dipole_vector':mu,'forces':forces}
        return results


    def _kernel_forces(self,mol,q,deriv="numerical"):
        dKdr = self.Kernel._an_kernel_gradient(mol)
        #dKdr = self.Kernel._num_kernel_gradient(mol)
        nAt  = len(q)
        f_k  = np.zeros((nAt,3))
        for j in range(nAt):
            for direction in range(3):
                eneg_drj = np.matmul(dKdr[j][direction],self.weights)
                for i in range(nAt):
                    f_k[j,direction] += -q[i]*eneg_drj[i]
        return f_k
    
    
    def _build_A(self,mols):
        all_as = []
        dim = 0
        for mol in mols:
            qe   = charge_eq(mol,scale_atsize=self.scale_atsize,radius_type=self.radius_type, hard=self.hard_lib)
            dim += qe.nAt
            all_as.append(qe.get_A())
        A = _block_diag(all_as,dim)
        return A


    def _build_A_X_O(self,mols=None,systems_charges=None,nAtoms = None):
        """Computes blocked matrixes needed for predict function

        Parameters
        ----------
            predict_set : list
                List of atoms objects
            predict_system_charges : list
                List of integers representing charges of system
            nAtoms : int
                Number of atoms in the prediction set
        
        Returns
        -------
            A : array
                Blocked inverse hardness matrix
            X : array
                Transformation matrix for N+1 vector 
            O : array
                Vector of zeroes and -Q of the system
            R : array
                Blocked dipole transformation matrix
            Xback : array
                Transformation matrix for creating N vector
        """   

        if nAtoms == None:
            nAtoms = 0
            for mol in mols:
                nAtoms += len(mol)
        count_row = 0
        count_col = 0
        X = np.zeros((nAtoms + len(mols),nAtoms))
        O = np.zeros((nAtoms + len(mols)))
        countChar = 0
        dim = 0
        all_as = []
        for mol in mols:
            qe   = charge_eq(mol,scale_atsize=self.scale_atsize,radius_type=self.radius_type, hard=self.hard_lib)
            dim += qe.nAt+1
            all_as.append(qe.get_Abar())
            for at in mol:
                X[count_row,count_col] = 1
                count_row += 1
                count_col += 1
            O[count_row] = -systems_charges[countChar]
            countChar += 1
            count_row +=1
        A = _block_diag(all_as,dim)
        return A, X, O



    def _build_prediction(self,mols,systems_charges,nAtoms = None):
        """Computes blocked matrixes needed for predict function

        Parameters
        ----------
            predict_set : list
                List of atoms objects
            predict_system_charges : list
                List of integers representing charges of system
            nAtoms : int
                Number of atoms in the prediction set
        
        Returns
        -------
            A : array
                Blocked inverse hardness matrix
            X : array
                Transformation matrix for N+1 vector 
            O : array
                Vector of zeroes and -Q of the system
            R : array
                Blocked dipole transformation matrix
            Xback : array
                Transformation matrix for creating N vector
        """   

        all_rs = []
        if nAtoms == None:
            nAtoms = 0
            for mol in mols:
                nAtoms += len(mol)
        count_row = 0
        count_col = 0
        count_rowBack = 0
        count_colBack = 0
        X = np.zeros((nAtoms + len(mols),nAtoms))
        Xback = np.zeros((nAtoms,nAtoms+len(mols)))
        O = np.zeros((nAtoms + len(mols)))
        countChar = 0
        dim = 0
        all_as = []
        for mol in mols:
            all_rs.append(_get_R(mol))
            qe   = charge_eq(mol,scale_atsize=self.scale_atsize,radius_type=self.radius_type, hard=self.hard_lib)
            dim += qe.nAt+1
            all_as.append(qe.get_Abar())
            for at in mol:
                Xback[count_rowBack,count_colBack] = 1
                X[count_row,count_col] = 1
                count_row += 1
                count_col += 1
                count_colBack += 1
                count_rowBack += 1
            O[count_row] = -systems_charges[countChar]
            countChar += 1
            count_row +=1
            count_colBack += 1
        A = _block_diag(all_as,dim)
        R = _block_diag_rect(all_rs,dim,len(all_rs)*3)
        return A, X, O, R, Xback

    def calculateEnergy(self, mol, charge=0):
        """Calculates charges, energy, forces and dipole_vector for single atoms object.

        Parameters
        ----------
            mol : obj
                Atoms object
            charge : int
                Charge of mol

        Returns
        -------
            results: dict
                Dictionary of results

        """

        K = self.Kernel.kernel_matrix(mol_set1 = [mol],kerneltype="predicting")#,predict_set=self.training_set,kerneltype="predicting") 
        eneg   = np.matmul(K,self.weights)
        qe     = charge_eq(mol,Qtot=charge,eneg=eneg,scale_atsize=self.scale_atsize,DoLJ=self.LJ,radius_type=self.radius_type, hard=self.hard_lib)
        qe.calc_Charges()
        qe.calc_Eqeq()
        energy  = qe.E 
        return energy*Hartree
    

    def save_kQEq(self):
        np.save("weights.npy", self.weights)
        np.save("sparseSOAP.npy", self.Kernel.representing_set_descriptors)

        #np.savetxt("weights.kqeq", self.weights, delimiter="\n")
        #np.savetxt("sparseSOAP.kqeq", self.Kernel.representing_set_descriptors, delimiter="\n")
    
    def load_kQEq(self):
        if self.sparse == True:
            self.weights = np.load("weights.npy")
            self.Kernel.representing_set_descriptors = np.load("sparseSOAP.npy")
        elif self.sparse == False:
            self.weights = np.load("weights.npy")
            self.Kernel.training_set_descriptors = np.load("sparseSOAP.npy")
