from pickletools import optimize
import numpy as np
from ase.io import *
from ase.data import covalent_radii
from ase.units import Bohr
import random
from quippy.descriptors import Descriptor
from numba import njit
from dscribe.descriptors import SOAP



class SOAPKernel():
    def __init__(self,multi_SOAP=False,descriptor_dict=None,
                training_set=None, training_system_charges = None,
                validation_set=None,validation_system_charges = None):
        """
        Parameters
        ----------
        Kernel : string
            Kernel function
        Descriptor : string
            descriptor
        descriptor_dict : dict
            dictionary of descriptor specific hyperparameters
        training_set : list
            List of atoms objects containing training data
        training_system_charges : list
            List of charges of training data 
        validation_set : list 
            (Optional) List of atoms objects containing validation data. Used to precompute matrices to accellerate hyperparameter optimization.
        validation_system_charges : list
            List of charges of validation data
        """
        self.sparse=False
        self.descriptor_dict = descriptor_dict
        self.multi_SOAP = multi_SOAP
        if training_set is not None:
            self.training_set = training_set
            self.training_set_descriptors = self.descriptor_atoms(training_set)
            if training_system_charges == None:
                self.training_system_charges = [0 for temp in training_set]   
            else: 
                self.training_system_charges = training_system_charges
        else:
            self.training_set_descriptors = None
            self.training_set = None
            self.training_system_charges = None
        if validation_set is not None:
            self.validation_set_descriptors = self.descriptor_atoms(validation_set)
        else:
            self.validation_set_descriptors = None

        self.validation_set = validation_set
        if validation_system_charges == None and validation_set is not None:
            self.validation_system_charges = [0 for temp in validation_set]
        else: 
            self.validation_system_charges = validation_system_charges

    def _create_representationFPS(self, sparse_positions):
        """
        Function for creation of representative set based on FPS
        
        
        Parameters
        ----------
        sparse_positions : list
            data picked by FPS from kernel_qeq class
        """
        if self.multi_SOAP == False: 
            self.representing_set_descriptors = [self.training_set_descriptors[pos] for pos in sparse_positions]
            self.sparse=True
        elif self.multi_SOAP == True:
            self.representing_set_descriptors = []
            for soap_n,sopa_d in enumerate(self.training_set_descriptors):
                rep_one = [self.training_set_descriptors[soap_n][pos] for pos in sparse_positions]
                self.representing_set_descriptors.append(rep_one)
            self.sparse=True

    def _create_representationCUR(self,sparse_count):
        """
        Function for selection and creation of representative set based on CUR decomposition 
        
        
        Parameters
        ----------
        sparse_count : 
            number of data points in represetative set
            
        Returns
        -------
        picked : list 
            list of index for representative set 
            
        """
        np.random.seed(10)
        k = 10
        epsilon = 0.45
        c = k*(np.log(k))/(epsilon**2)
        while sparse_count > c:
            k += 3
            c = k*(np.log(k))/(epsilon**2)
        if self.multi_SOAP == False:
            n_struc = len(self.training_set_descriptors)
            trainT = np.transpose(self.training_set_descriptors)
        elif self.multi_SOAP == True:
            n_struc = len(self.training_set_descriptors[0])
            trainT = np.transpose(self.training_set_descriptors[0])
        u,s,V = np.linalg.svd(trainT)
        V = V[:k,:]
        probV = 1/k*np.sum(V**2,axis=0)
        C_vec = [np.random.choice(2, p=[1 - min(1, c * probV[strucN]), min(1, c * probV[strucN])],replace=False) for strucN in range(n_struc)]
        picked = np.nonzero(C_vec)[0]
        if self.multi_SOAP == False:
            self.representing_set_descriptors = [self.training_set_descriptors[pos] for pos in picked]
        elif self.multi_SOAP == True:
            self.representing_set_descriptors = []
            for ind in range(len(self.training_set_descriptors)):
                self.representing_set_descriptors.append([self.training_set_descriptors[ind][pos] for pos in picked])
        self.sparse = True
        return picked
        
    def descriptor_atoms(self,molecules):
        """
        Yields descriptors for the input molecules based on the choosen descriptor.
        Descriptor uses the hyperparameter dictionary.
        
        
        Parameters
        ----------
        molecules : list
            list of ase atoms objects
            
        Returns
        -------
        descriptor : list 
            elements of the list depend on the choosen descriptor
            
        """ 
        if self.multi_SOAP:
            nmax =  self.descriptor_dict['nmax']
            lmax = self.descriptor_dict['lmax']
            rcut = self.descriptor_dict['rcut']
            if self.training_set is not None:
                species = np.unique([spec  for mol in self.training_set for spec in  mol.get_chemical_symbols()])
            else:
                species = np.unique([spec  for mol in molecules for spec in  mol.get_chemical_symbols()])
            periodic = self.descriptor_dict['periodic']
            sigma = self.descriptor_dict['sigma']
            descriptor = [self.atomic_descriptor_SOAP(molecules,species,rcut[soap_ind],lmax[soap_ind],nmax[soap_ind],periodic=periodic[soap_ind],sigma=sigma[soap_ind]) for soap_ind in range(len(nmax))]
        else:
            nmax =  self.descriptor_dict['nmax']
            lmax = self.descriptor_dict['lmax']
            rcut = self.descriptor_dict['rcut']
            if self.training_set is not None:
                species = np.unique([spec  for mol in self.training_set for spec in  mol.get_chemical_symbols()])
            else:
                species = np.unique([spec  for mol in molecules for spec in  mol.get_chemical_symbols()])
            periodic = self.descriptor_dict['periodic']
            sigma = self.descriptor_dict['sigma']
            descriptor =  self.atomic_descriptor_SOAP(molecules,species,rcut,lmax,nmax,periodic=False,sigma=sigma)
        return descriptor
            
            

    
    def atomic_descriptor_SOAP(self,molecules,species,rcut,lmax,nmax,periodic=False,sigma=1.0):
        """
        SOAP encodes regions of atomic geometries by using a local expansion of a gaussian smeared atomic density with orthonormal functions based on spherical harmonics and radial basis functions.
        
        Parameters
        ----------
        molecules :  list
            list of ase atoms objects
        species : list
             The chemical species as a list of atomic numbers or as a list of chemical symbols.
        rcut : float
            A cutoff for local region in angstroms. Should be bigger than 1 angstrom.
        nmax : int  
            The number of radial basis functions.
        lmax : int 
            The maximum degree of spherical harmonics.
        periodic : bool 
            Determines whether the system is considered to be periodic
        sigma : float
            The standard deviation of the gaussians used to expand the atomic density.
            
        Returns
        -------
        soap_atoms : list
            list of normalized soap vectors.
        """
        
        soap = SOAP(
            species=species,
            periodic=periodic,
            rcut=rcut,  
            nmax=nmax,
            sigma=sigma,
            lmax=lmax)
        soap_atoms=[]
        for mol in molecules:
            soap_vec = soap.create(mol)
            normal_soap_vector = []
            for sv in soap_vec:
                
                norm = np.linalg.norm(sv)
                normal_soap_vector.append(sv/norm)
            soap_atoms.extend(normal_soap_vector) 
        return soap_atoms    


    
    def soap_kernel(self,mol_set1,mol_set2,zeta=2,kerneltype='general'):
        """
        The SOAP kernel between two atomic environments can be retrieved as a normalized polynomial kernel of the partial powers         spectrums.
        
        Parameters
        ----------
        molsA : list
            list of ase atoms objects
        zeta : float
            exponent
            
        Returns
        -------
        K : matrix
            Kernel Matrix. Matrix elements range from 0 to 1
        """
        if self.multi_SOAP:
            deltas = self.descriptor_dict['deltas']
            if kerneltype=='general':
                if mol_set2 == None:
                    desc1 = self.descriptor_atoms(mol_set1)
                    K_multi = []
                    for soap_ind,soap_des in enumerate(desc1):
                        k = np.matmul(soap_des,np.transpose(soap_des))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                    K = 0
                    for i in K_multi:
                        K=np.add(i,K)        
                        k = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                    K = 0
                    for i in K_multi:
                        K=np.add(i,K)
                return K
            elif kerneltype=='training':
                if self.sparse == False:
                    desc1 = self.training_set_descriptors
                    K_multi = []
                    for soap_ind,soap_des in enumerate(desc1):
                        k = np.matmul(soap_des,np.transpose(soap_des))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                    K = 0
                    for i in K_multi:
                        K=np.add(i,K)
                    return K
                elif self.sparse == True:
                    K_NM_multi = []
                    K_MM_multi = []
                    desc1 = self.training_set_descriptors
                    desc2 = self.representing_set_descriptors
                    for soap_ind,soap_des in enumerate(desc1): 
                        k_nm = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_NM_multi.append(np.multiply(deltas[soap_ind],k_nm))
                        k_mm = np.matmul(desc2[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_MM_multi.append(np.multiply(deltas[soap_ind],k_mm))
                    K_NM = 0
                    K_MM = 0
                    for i in K_NM_multi:
                        K_NM=np.add(i,K_NM)
                    for i in K_MM_multi:
                        K_MM=np.add(i,K_MM)
                    return K_NM,K_MM
            elif kerneltype=='validation':
                desc1 = self.validation_set_descriptors
                # desc2 = self.training_set_descriptors
                K_multi = []
                if self.sparse == False:
                    desc2 = self.training_set_descriptors
                    for soap_ind,soap_des in enumerate(desc1): 
                        k = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                if self.sparse == True:
                    desc2 = self.representing_set_descriptors
                    for soap_ind,soap_des in enumerate(desc1): 
                        k = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                #for soap_ind,soap_des in enumerate(desc1):
                #    k = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                #    K_multi.append(np.multiply(deltas[soap_ind],k))
                K = 0
                for i in K_multi:
                    K=np.add(i,K)
                return K 
            elif kerneltype=='predicting':
                desc1 = self.descriptor_atoms(mol_set1)
                K_multi = []
                if self.sparse == False:
                    
                    desc2 = self.training_set_descriptors
                    for soap_ind,soap_des in enumerate(desc1): 
                        k = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                    K = 0
                    for i in K_multi:
                        K=np.add(i,K)
                    return K
                if self.sparse == True:   
                    desc2 = self.representing_set_descriptors
                    for soap_ind,soap_des in enumerate(desc1): 
                        k = np.matmul(desc1[soap_ind],np.transpose(desc2[soap_ind]))**zeta
                        K_multi.append(np.multiply(deltas[soap_ind],k))
                    K = 0
                    for i in K_multi:
                        K=np.add(i,K)
                    return K
            
        else:
            if kerneltype=='general':
                if mol_set2 == None:
                    desc1 = self.descriptor_atoms(mol_set1)
                    dim1  = len(desc1)
                    K = np.matmul(desc1,np.transpose(desc1))**zeta
                else:
                    desc1 = self.descriptor_atoms(mol_set1)
                    desc2 = self.descriptor_atoms(mol_set2)
                    dim1  = len(desc1)
                    dim2  = len(desc2)
                    K = np.matmul(desc1,np.transpose(desc2))**zeta
                return K
            elif kerneltype=='predicting':
                desc1 = self.descriptor_atoms(mol_set1)
                if self.sparse == False:
                    K = np.matmul(desc1,np.transpose(self.training_set_descriptors))**zeta
                if self.sparse == True:
                    K = np.matmul(desc1,np.transpose(self.representing_set_descriptors))**zeta
                return K
            elif kerneltype=='training':
                if self.sparse == False:
                    K = np.matmul(self.training_set_descriptors,np.transpose(self.training_set_descriptors))**zeta
                    return K
                elif self.sparse == True:
                    K_nm = np.matmul(self.training_set_descriptors,np.transpose(self.representing_set_descriptors))**zeta
                    K_mm = np.matmul(self.representing_set_descriptors,np.transpose(self.representing_set_descriptors))**zeta
                    return K_nm, K_mm

            elif kerneltype=='validation':
                if self.sparse == False:
                    K = np.matmul(self.validation_set_descriptors,np.transpose(self.training_set_descriptors))**zeta
                if self.sparse == True:
                    K = np.matmul(self.validation_set_descriptors,np.transpose(self.representing_set_descriptors))**zeta
                return K
        
    
    
    def kernel_matrix(self,mol_set1=None,mol_set2=None,kerneltype='general'):
        """
        Returns the Kernel matrix for the specified kernel.
        
        Parameters
        ----------
        training_set : list
            list of ase atoms objects
        predict_set: list
            list of ase atoms objects
        
        Returns
        -------
        K : matrix
            Kernel Matrix
        """

        K = self.soap_kernel(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
        return K
    
    
    def calculate_function(self,mol_set1,zeta=2):
        desc1,der1 = self._descr_deriv_atoms([mol_set1],deriv="numerical")
        desc_train = self.representing_set_descriptors
        normal_soap_vector = []
        norms = []
        for sv in desc1:
            norm = np.linalg.norm(sv)
            normal_soap_vector.append(sv/norm)
            norms.append(norm)
        norms = np.array(norms)
        norms2 = norms**2
        K_final = np.matmul(normal_soap_vector,np.transpose(desc_train))**zeta
        temp1 = []
        temp2 = []
        deriv_final = []
        for a in range(len(der1)):
            for b in range(3):
                for c in range(len(der1)):
                    temp1.append(der1[c][a][b])
                    #temp1.append(der1[a][c][b])
                temp2.append(temp1)
                temp1 = []
            deriv_final.append(temp2) 
            temp2 = []
        dKdr = [] 
        deriv_final = np.array(deriv_final)
        hold = (np.einsum('ij,arkj->arik',desc1,deriv_final,optimize="greedy"))
        vd = np.einsum('akii->aki',hold,optimize="greedy")
        r1 = np.einsum('akl,lj,l->aklj',vd,desc1,1/norms,optimize="greedy")#"greedy")
        f1 = np.einsum('aklj,l->aklj',deriv_final,norms,optimize="greedy")
        normDer = np.einsum('aklj,l->aklj',f1-r1,1/norms2,optimize="greedy")
        #t1 = np.einsum('ij,lj->il',normal_soap_vector,desc_train,optimize="greedy")**(zeta-1)
        #t2 = np.einsum('akij,lj->akil',normDer,desc_train,optimize="greedy")
        #dKdr = zeta*np.einsum('akij,ij->akij',t2,t1,optimize="greedy")
#        t1 = np.matmul(normal_soap_vector,np.transpose(desc_train))**(zeta-1)
        t1 = np.einsum('ij,lj->il',normal_soap_vector,desc_train,optimize="greedy")**(zeta-1)
        dKdr = zeta*np.einsum('akij,lj,il->akil',normDer,desc_train,t1,optimize="greedy")        
        return K_final,dKdr
    
    
    def _soap_kernel_descriptors(self,ref_mol,descriptors_training,zeta=2):
        """
        The SOAP kernel between two atoms object and precomputed descriptors (for gradient code)
        
        Parameters
        ----------
        molA : obj
            ase atoms objects
        zeta : float
            exponent
        descriptors_training : list
            precomputed SOAP descriptors
            
        Returns
        -------
        K : matrix
            Kernel Matrix. Matrix elements range from 0 to 1
        """
        if self.multi_SOAP:
            deltas = self.descriptor_dict['deltas']
            desc1 = self.descriptor_atoms(ref_mol)
            dim1  = len(desc1)
            dim2  = len(descriptors_training)
            K_multi = []
            for soap_ind,soap_des in enumerate(desc1): 
                k = np.matmul(desc1[soap_ind],np.transpose(descriptors_training[soap_ind]))**zeta
                K_multi.append(np.multiply(deltas[soap_ind],k))
            K = 0
            for i in K_multi:
                K=np.add(i,K)
            return K
        else:
            desc1 = self.descriptor_atoms(ref_mol)
            #desc2 = self.descriptor_atoms(compare_set)
            #dim1  = len(desc1)
            dim2  = len(descriptors_training)
            K = np.matmul(desc1,np.transpose(descriptors_training))**zeta
            return K

    def _an_kernel_gradient(self,ref_mol,zeta=2):
        if self.sparse == False:
#                desc_train = self.descriptor_atoms(training_set)
            desc_train = self.training_set_descriptors
        elif self.sparse == True:
            desc_train = self.representing_set_descriptors
        desc1,der1 = self._descr_deriv_atoms([ref_mol],deriv="numerical")
        temp1 = []
        temp2 = []
        deriv_final = []
        if self.multi_SOAP == False:
        # loop for reorganizing dscribed ordering to original numerical one
            for a in range(len(der1)):
                for b in range(3):
                    for c in range(len(der1)):
                        temp1.append(der1[c][a][b])
                        #temp1.append(der1[a][c][b])
                    temp2.append(temp1)
                    temp1 = []
                deriv_final.append(temp2) 
                temp2 = []
                normal_soap_vector = []
                norms = []
                for sv in desc1:
                    norm = np.linalg.norm(sv)
                    normal_soap_vector.append(sv/norm)
                    norms.append(norm)
                norms = np.array(norms)
            dKdr = [] 
            for iatom in range(len(ref_mol)):
                dKdri = []
                for direction in range(3):
                    inter = deriv_final[iatom][direction] # derivation of atom in 3 direction based on all other atoms
                    tt1 = (np.matmul(desc1,np.transpose(inter))) 
                    vd = np.diag(tt1)
                    tt2 = np.multiply(desc1,vd[:,None])
                    r1 = tt2/norms[:,None] # each atomistic SOAP descriptor is normalized 
                    f1 = np.multiply(inter,norms[:,None]) 
                    norms2 = norms**2
                    normDer=(f1-r1)/norms2[:,None]
                    t1 = (np.matmul(normal_soap_vector,np.transpose(desc_train))**(zeta-1))
                    t2 = (np.matmul(normDer,np.transpose(desc_train)))
                    K = zeta*np.multiply(t1,t2)
                    dKdri.append(K)
                dKdr.append(dKdri)
            return dKdr
        elif self.multi_SOAP == True:
            deltas = self.descriptor_dict['deltas']
            norms_final = []
            normal_soap_vector = []
            for tempMAIN in der1:
                deriv_temp = []
                for a in range(len(tempMAIN)):
                    for b in range(3):
                        for c in range(len(tempMAIN)):
                            temp1.append(tempMAIN[c][a][b])
                            #temp1.append(der1[a][c][b])
                        temp2.append(temp1)
                        temp1 = []
                    deriv_temp.append(temp2) 
                    temp2 = []
                    
                    
                    normsAll = []
                    for descTemp in desc1:
                        norms = []
                        normal_soap_vector_temp = []
                        for sv in descTemp:
                            norm = np.linalg.norm(sv)
                            normal_soap_vector_temp.append(sv/norm)
                            norms.append(norm)
                        normal_soap_vector.append(normal_soap_vector_temp)
                        norms_final.append(np.array(norms))
                deriv_final.append(deriv_temp)
            dKtemp = []
            
            for soap_ind,soap_des in enumerate(desc1):
                dKdr = [] 
                for iatom in range(len(ref_mol)):
                    dKdri = []
                    for direction in range(3):
                        inter = deriv_final[soap_ind][iatom][direction]
                        tt1 = (np.matmul(desc1[soap_ind],np.transpose(inter))) 
                        vd = np.diag(tt1)
                        tt2 = np.multiply(desc1[soap_ind],vd[:,None])
                        nm_b = norms_final[soap_ind][:,None]
                        r1 = tt2/nm_b # each atomistic SOAP descriptor is normalized 
                        f1 = np.multiply(inter,norms_final[soap_ind][:,None]) 
                        norms2 = norms_final[soap_ind]**2
                        normDer=(f1-r1)/norms2[:,None]                                       
                        t1 = (np.matmul(normal_soap_vector[soap_ind],np.transpose(desc_train[soap_ind]))**(zeta-1))
                        t2 = (np.matmul(normDer,np.transpose(desc_train[soap_ind])))
                        k = zeta*np.multiply(t1,t2)
                        K = np.multiply(deltas[soap_ind],k)
                        dKdri.append(K)
                    dKdr.append(dKdri)
                dKtemp.append(dKdr)
            dKf = 0
            for i in dKtemp:
                dKf=np.add(i,dKf)
            return dKf


   

    def _descr_deriv_atoms(self,molecules,deriv):
        nmax =  self.descriptor_dict['nmax']
        lmax = self.descriptor_dict['lmax']
        rcut = self.descriptor_dict['rcut']
        if self.training_set is not None:
            species = np.unique([spec  for mol in self.training_set for spec in  mol.get_chemical_symbols()])
        else:
            species = np.unique([spec  for mol in molecules for spec in  mol.get_chemical_symbols()])
        periodic = self.descriptor_dict['periodic']
        sigma = self.descriptor_dict['sigma']
        if self.multi_SOAP == False:
            descriptor, derivative =  self._atomic_desc_deriv_SOAP(molecules,species,rcut,lmax,nmax,deriv,periodic=False,sigma=sigma)
        elif self.multi_SOAP == True:
            derivative = []
            descriptor = []
            for soap_ind in range(len(nmax)):
                desctemp, dertemp  = self._atomic_desc_deriv_SOAP(molecules,species,rcut[soap_ind],lmax[soap_ind],nmax[soap_ind],deriv,periodic=periodic[soap_ind],sigma=sigma[soap_ind])
                derivative.append(dertemp)
                descriptor.append(desctemp)
                
        return np.array(descriptor), np.array(derivative)

    def _atomic_desc_deriv_SOAP(self,molecule,species,rcut,lmax,nmax,deriv,periodic=False,sigma=1.0):
        soap = SOAP(
            species=species,
            periodic=periodic,
            rcut=rcut,
            nmax=nmax,
            sigma=sigma,
            lmax=lmax) 
        soap_d, soap_vec = soap.derivatives(molecule,method="numerical",attach=True) #analytical derivatives in SOAPs are not usable with 1.2.0 version
        return soap_vec, soap_d   



class SOAPKernelQUIP():
    def __init__(self,Kernel='Delta',Descriptor='Delta',
                multi_SOAP=False,descriptor_dict=None,
                training_set=None, training_system_charges = None,
                validation_set=None,validation_system_charges = None):
        self.Kernel = Kernel
        self.sparse=False
        self.Descriptor = Descriptor
        #self.descriptor_dict = descriptor_dict
        self.multi_SOAP = multi_SOAP
        self.representing_set_descriptors = None
        #self.Descriptor_description = f"soap cutoff={self.descriptor_dict['rcut']} zeta=2 delta = 1.0 atom_sigma = {self.descriptor_dict['sigma']} l_max={self.descriptor_dict['lmax']} n_max={self.descriptor_dict['nmax']} add_species=F"
        self.Descriptor_description = descriptor_dict
        if training_set is not None:
            self.training_set = training_set
            self.training_set_descriptors = self.atomic_descriptor_SOAP(training_set)
        else:
            self.training_set_descriptors = None
            self.training_set = None
        if validation_set is not None:
            self.validation_set_descriptors = self.atomic_descriptor_SOAP(validation_set)
        else:
            self.validation_set_descriptors = None
        if training_system_charges == None:
            self.training_system_charges = [0 for temp in training_set]   
        else: 
            self.training_system_charges = training_system_charges
        self.validation_set = validation_set
        if validation_system_charges == None and validation_set is not None:
            self.validation_system_charges = [0 for temp in validation_set]
        else: 
            self.validation_system_charges = validation_system_charges
        
    def atomic_descriptor_SOAP(self,molecules):
        #self.Descriptor_description = f"soap cutoff={rcut} zeta=2 delta = 1.0 atom_sigma = {sigma} l_max={lmax} n_max={nmax} add_species=F"
        soap = Descriptor(self.Descriptor_description)
        soaps = soap.calc(molecules)
        soap_atoms=[]
        for mol in soaps:
            for at in mol["data"]:
                #norm = np.linalg.norm(at)
                soap_atoms.append(at)
        return soap_atoms

 
    def kernel_matrix(self,mol_set1=None,mol_set2=None,kerneltype='general'):
        K = self.soap_kernel(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
        return K

    
    def soap_kernel(self,mol_set1,mol_set2,zeta=2,kerneltype='general'):
        if kerneltype=='general':
            if mol_set2 == None:
                desc1 = self.atomic_descriptor_SOAP(mol_set1)
                dim1  = len(desc1)
                K = np.matmul(desc1,np.transpose(desc1))**zeta
            else:
                desc1 = self.atomic_descriptor_SOAP(mol_set1)
                desc2 = self.atomic_descriptor_SOAP(mol_set2)
                dim1  = len(desc1)
                dim2  = len(desc2)
                K = np.matmul(desc1,np.transpose(desc2))**zeta
            return K
        elif kerneltype=='predicting':
            desc1 = self.atomic_descriptor_SOAP(mol_set1)
            if self.sparse == False:
                K = np.matmul(desc1,np.transpose(self.training_set_descriptors))**zeta
            if self.sparse == True:
                K = np.matmul(desc1,np.transpose(self.representing_set_descriptors))**zeta
            return K
        elif kerneltype=='training':
            if self.sparse == False:
                K = np.matmul(self.training_set_descriptors,np.transpose(self.training_set_descriptors))**zeta
                return K
            elif self.sparse == True:
                K_nm = np.matmul(self.training_set_descriptors,np.transpose(self.representing_set_descriptors))**zeta
                K_mm = np.matmul(self.representing_set_descriptors,np.transpose(self.representing_set_descriptors))**zeta
                return K_nm, K_mm
    
    def _create_representationCUR(self,sparse_count):
        """
        Function for selection and creation of representative set based on CUR decomposition 
        
        
        Parameters
        ----------
        sparse_count : int
            number of data points in represetative set
            
        Returns
        -------
        picked : list 
            list of index for representative set 
            
        """
        np.random.seed(10)
        k = 10
        epsilon = 0.45
        c = k*(np.log(k))/(epsilon**2)
        while sparse_count > c:
            k += 3
            c = k*(np.log(k))/(epsilon**2)
        if self.multi_SOAP == False:
            n_struc = len(self.training_set_descriptors)
            trainT = np.transpose(self.training_set_descriptors)
        elif self.multi_SOAP == True:
            n_struc = len(self.training_set_descriptors[0])
            trainT = np.transpose(self.training_set_descriptors[0])
        u,s,V = np.linalg.svd(trainT)
        V = V[:k,:]
        probV = 1/k*np.sum(V**2,axis=0)
        picked = np.argsort(probV)[:sparse_count]
        if self.multi_SOAP == False:
            self.representing_set_descriptors = [self.training_set_descriptors[pos] for pos in picked]
        elif self.multi_SOAP == True:
            self.representing_set_descriptors = []
            for ind in range(len(self.training_set_descriptors)):
                self.representing_set_descriptors.append([self.training_set_descriptors[ind][pos] for pos in picked])
        self.sparse = True
        return picked


    def _descr_deriv_atoms(self,molecules):
        soap = Descriptor(self.Descriptor_description)
        soaps = soap.calc(molecules, grad=True)
        soap_atoms=[]
        soap_derivatives = []
        neigh_id = []
        for mol in soaps:
            for at in mol["data"]:
                #norm = np.linalg.norm(at)
                soap_atoms.append(at)
            for at in mol["grad_data"]:
                soap_derivatives.append(at)
            for at in mol["grad_index_0based"]:    
                neigh_id.append(at)        
        return soap_atoms, soap_derivatives, neigh_id

    def _num_kernel_gradient(self,ref_mol,h=0.001):
        if self.sparse == False:
            desc_train = self.training_set_descriptors
        elif self.sparse == True:
            desc_train = self.representing_set_descriptors
        dKdr = []
        for i in range(len(ref_mol)):
            dKdri = []
            for direction in range(3):
                dKdri.append(self._num_gradient(ref_mol,h=h,direction=direction,iatom=i))
            dKdr.append(dKdri)
        return dKdr
    
    def _an_kernel_gradient(self,ref_mol,zeta=2):
        if self.sparse == False:
    #                desc_train = self.descriptor_atoms(training_set)
            desc_train = self.training_set_descriptors
        elif self.sparse == True:
            desc_train = self.representing_set_descriptors
        desc1,der1, id = self._descr_deriv_atoms([ref_mol])
        id = [list(item) for item in id]
        derAll = []
        derOne = [[],[],[]]
        lenDer = len(der1[0][0])
        lib_der = {}
        for a in range(len(der1)):
            lib_der[str(id[a])] = der1[a] 
        for a1 in range(len(ref_mol)):
            for a2 in range(len(ref_mol)):
                if [a1,a2] in id:
                    derOne[0].append(lib_der[str([a1,a2])][0])
                    derOne[1].append(lib_der[str([a1,a2])][1])
                    derOne[2].append(lib_der[str([a1,a2])][2])
                else:
                    for di in range(3):
                        derOne[di].append(np.zeros(lenDer))
            derAll.append(derOne)
            derOne = [[],[],[]]
        dKdr = []
        temp1 = []
        temp2 = []
        deriv_final = []
        for a in range(len(derAll)):
            for b in range(3):
                for c in range(len(derAll)):
                    temp1.append(derAll[c][b][a])
                temp2.append(temp1)
                temp1 = []
            deriv_final.append(temp2) 
            temp2 = []
        for iatom in range(len(ref_mol)):
            dKdri = []
            for direction in range(3):
                curDer = (deriv_final[iatom][direction])
                t1 = (np.matmul(desc1,np.transpose(desc_train))**(zeta-1))
                t2 = np.matmul(curDer,np.transpose(desc_train))
                K = (zeta*np.multiply(t1,t2))
                dKdri.append(K)
            dKdr.append(dKdri)
        return dKdr

    def _num_gradient(self,mol,h=0.001,direction=0,iatom=0):
        tmpmol = mol.copy()
        pos = tmpmol.get_positions()
        #pos /= Bohr
        pos[iatom][direction] += h
        tmpmol.set_positions(pos)
        Kplus = self._soap_kernel_descriptors([tmpmol])
        pos[iatom][direction] += -2.0*h
        tmpmol.set_positions(pos)
        Kminus = self._soap_kernel_descriptors([tmpmol])
        pos[iatom][direction] += h
        tmpmol.set_positions(pos)

        return (Kplus-Kminus)/(2.0*h)


    def _soap_kernel_descriptors(self,ref_mol,zeta=2):
        if self.sparse == False:
            desc_train = self.training_set_descriptors
        elif self.sparse == True:
            desc_train = self.representing_set_descriptors
        desc1 = self.atomic_descriptor_SOAP(ref_mol)
        K = np.matmul(desc1,np.transpose(desc_train))**zeta
        return K

class soapGAPkQEqQUIP():
    
    def __init__(self,multi_SOAP=False,descriptor_dict_GAP=None, descriptor_dict_kQEq=None,
                training_set=None, training_system_charges = None,
                validation_set=None,validation_system_charges = None):
        self.sparse=False
        self.multi_SOAP = multi_SOAP
        self.representing_set_descriptors_GAP = None
        self.representing_set_descriptors_kQEq = None
        self.Descriptor_description_GAP = descriptor_dict_GAP
        self.Descriptor_description_kQEq = descriptor_dict_kQEq

        if training_set is not None:
            self.training_set = training_set
            self.training_set_descriptors_GAP = self.atomic_descriptor_SOAP(training_set, self.Descriptor_description_GAP)
            self.training_set_descriptors_kQEq =  self.atomic_descriptor_SOAP(training_set, self.Descriptor_description_kQEq)
        else:
            self.training_set_descriptors_GAP = None
            self.training_set_descriptors_kQEq = None
            self.training_set = None
        if validation_set is not None:
            self.validation_set = validation_set
            self.validation_set_descriptors_GAP = self.atomic_descriptor_SOAP(validation_set, self.Descriptor_description_GAP)
            self.validation_set_descriptors_kQEq =  self.atomic_descriptor_SOAP(validation_set, self.Descriptor_description_kQEq)
        else:
            self.validation_set_descriptors_GAP = None
            self.validation_set_descriptors_kQEq = None
        if training_system_charges == None:
            self.training_system_charges = [0 for temp in training_set]   
        else: 
            self.training_system_charges = training_system_charges
        self.validation_set = validation_set
        if validation_system_charges == None and validation_set is not None:
            self.validation_system_charges = [0 for temp in validation_set]
        else: 
            self.validation_system_charges = validation_system_charges
        
    def atomic_descriptor_SOAP(self,molecules,descript):
        #self.Descriptor_description = f"soap cutoff={rcut} zeta=2 delta = 1.0 atom_sigma = {sigma} l_max={lmax} n_max={nmax} add_species=F"
        soap = Descriptor(descript)
        soaps = soap.calc(molecules)
        soap_atoms=[]
        for mol in soaps:
            for at in mol["data"]:

                #norm = np.linalg.norm(at)
                soap_atoms.append(at)
        return soap_atoms


    def kernel_matrix_training(self,mol_set1=None,mol_set2=None):
        kerneltype='training'
        if self.sparse == True:
            K_nm_GAP,K_mm_GAP = self.soap_kernel_GAP(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
            K_nm_kQEq,K_mm_kQEq = self.soap_kernel_kQEq(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
            return K_nm_GAP,K_mm_GAP,K_nm_kQEq,K_mm_kQEq 
        elif self.sparse == False:
            K_GAP = self.soap_kernel_GAP(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
            K_kQEq = self.soap_kernel_kQEq(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
            return K_GAP,K_kQEq


    def kernel_matrix_prediction(self,mol_set1=None,mol_set2=None):
        kerneltype='predicting'
        K_GAP = self.soap_kernel_GAP(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
        K_kQEq = self.soap_kernel_kQEq(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
        return K_GAP,K_kQEq

    def soap_kernel_GAP(self,mol_set1,mol_set2,zeta=2,kerneltype='general'):

        if kerneltype=='general':
            if mol_set2 == None:
                desc1 = self.atomic_descriptor_SOAP(mol_set1,self.Descriptor_description_GAP)
                dim1  = len(desc1)
                K = np.matmul(desc1,np.transpose(desc1))**zeta
            else:
                desc1 = self.atomic_descriptor_SOAP(mol_set1,self.Descriptor_description_GAP)
                desc2 = self.atomic_descriptor_SOAP(mol_set2,self.Descriptor_description_GAP)
                dim1  = len(desc1)
                dim2  = len(desc2)
                K = np.matmul(desc1,np.transpose(desc2))**zeta
            return K
        elif kerneltype=='predicting':
            desc1 = self.atomic_descriptor_SOAP(mol_set1,self.Descriptor_description_GAP)
            if self.sparse == False:
                K = np.matmul(desc1,np.transpose(self.training_set_descriptors_GAP))**zeta
            if self.sparse == True:
                K = np.matmul(desc1,np.transpose(self.representing_set_descriptors_GAP))**zeta
            return K
        elif kerneltype=='training':
            if self.sparse == False:
                K = np.matmul(self.training_set_descriptors_GAP,np.transpose(self.training_set_descriptors_GAP))**zeta
                return K
            elif self.sparse == True:
                K_nm = np.matmul(self.training_set_descriptors_GAP,np.transpose(self.representing_set_descriptors_GAP))**zeta
                K_mm = np.matmul(self.representing_set_descriptors_GAP,np.transpose(self.representing_set_descriptors_GAP))**zeta
                return K_nm, K_mm
    
    def soap_kernel_kQEq(self,mol_set1,mol_set2,zeta=2,kerneltype='general'):

        if kerneltype=='general':
            if mol_set2 == None:
                desc1 = self.atomic_descriptor_SOAP(mol_set1,self.Descriptor_description_kQEq)
                dim1  = len(desc1)
                K = np.matmul(desc1,np.transpose(desc1))**zeta
            else:
                desc1 = self.atomic_descriptor_SOAP(mol_set1,self.Descriptor_description_kQEq)
                desc2 = self.atomic_descriptor_SOAP(mol_set2,self.Descriptor_description_kQEq)
                dim1  = len(desc1)
                dim2  = len(desc2)
                K = np.matmul(desc1,np.transpose(desc2))**zeta
            return K
        elif kerneltype=='predicting':
            desc1 = self.atomic_descriptor_SOAP(mol_set1,self.Descriptor_description_kQEq)
            if self.sparse == False:
                K = np.matmul(desc1,np.transpose(self.training_set_descriptors_kQEq))**zeta
            if self.sparse == True:
                K = np.matmul(desc1,np.transpose(self.representing_set_descriptors_kQEq))**zeta
            return K
        elif kerneltype=='training':
            if self.sparse == False:
                K = np.matmul(self.training_set_descriptors_kQEq,np.transpose(self.training_set_descriptors_kQEq))**zeta
                return K
            elif self.sparse == True:
                K_nm = np.matmul(self.training_set_descriptors_kQEq,np.transpose(self.representing_set_descriptors_kQEq))**zeta
                K_mm = np.matmul(self.representing_set_descriptors_kQEq,np.transpose(self.representing_set_descriptors_kQEq))**zeta
                return K_nm, K_mm


    
    def _create_representationCUR(self,sparse_count):
        """
        Function for selection and creation of representative set based on CUR decomposition 
        
        
        Parameters
        ----------
        sparse_count : int
            number of data points in represetative set
            
        Returns
        -------
        picked : list 9
            list of index for representative set 
            
        """
        np.random.seed(10)
        k = 10
        epsilon = 0.45
        c = k*(np.log(k))/(epsilon**2)
        while sparse_count > c:
            k += 3
            c = k*(np.log(k))/(epsilon**2)
        if self.multi_SOAP == False:
            n_struc = len(self.training_set_descriptors_GAP)
            trainT_GAP = np.transpose(self.training_set_descriptors_GAP)
            trainT_kQEq = np.transpose(self.training_set_descriptors_kQEq) 
        elif self.multi_SOAP == True:
            n_struc = len(self.training_set_descriptors_GAP[0])
            trainT_GAP = np.transpose(self.training_set_descriptors_GAP[0])
        u,s,V_gap = np.linalg.svd(trainT_GAP)
        u,s,V_kqeq = np.linalg.svd(trainT_kQEq)
        V_gap = V_gap[:k,:]
        V_kqeq = V_kqeq[:k,:]
        probV_gap = 1/k*np.sum(V_gap**2,axis=0)
        probV_kqeq = 1/k*np.sum(V_kqeq**2,axis=0)
        picked_gap = np.argsort(probV_gap)[:sparse_count]
        picked_kqeq = np.argsort(probV_kqeq)[:sparse_count]
        if self.multi_SOAP == False:
            self.representing_set_descriptors_GAP = [self.training_set_descriptors_GAP[pos] for pos in picked_gap]
            self.representing_set_descriptors_kQEq = [self.training_set_descriptors_kQEq[pos] for pos in picked_kqeq]
        elif self.multi_SOAP == True:
            self.representing_set_descriptors_GAP = []
            self.representing_set_descriptors_kQEq = []
            for ind in range(len(self.training_set_descriptors_GAP)):
                self.representing_set_descriptors_GAP.append([self.training_set_descriptors_GAP[ind][pos] for pos in picked_gap])
                self.representing_set_descriptors_kQEq.append([self.training_set_descriptors_kQEq[ind][pos] for pos in picked_kqeq])
        self.sparse = True
        return picked_gap, picked_kqeq


class DeltaKernel():
    def __init__(self, multi_SOAP=False,descriptor_dict=None,
                training_set=None, training_system_charges = None,
                validation_set=None,validation_system_charges = None):
        """
        Parameters
        ----------
        Kernel : string
            Kernel function
        Descriptor : string
            descriptor
        descriptor_dict : dict
            dictionary of descriptor specific hyperparameters
        training_set : list
            List of atoms objects containing training data
        training_system_charges : list
            List of charges of training data 
        validation_set : list 
            (Optional) List of atoms objects containing validation data. Used to precompute matrices to accellerate hyperparameter optimization.
        validation_system_charges : list
            List of charges of validation data
        """
        self.sparse=False
        self.Descriptor = Descriptor
        self.descriptor_dict = descriptor_dict
        self.multi_SOAP = multi_SOAP
        if training_set is not None:
            self.training_set = training_set
            self.training_set_descriptors = self.descriptor_atoms(training_set)
        else:
            self.training_set_descriptors = None
            self.training_set = None
        if validation_set is not None:
            self.validation_set_descriptors = self.descriptor_atoms(validation_set)
        else:
            self.validation_set_descriptors = None
        if training_system_charges == None:
            self.training_system_charges = [0 for temp in training_set]   
        else: 
            self.training_system_charges = training_system_charges
        self.validation_set = validation_set
        if validation_system_charges == None and validation_set is not None:
            self.validation_system_charges = [0 for temp in validation_set]
        else: 
            self.validation_system_charges = validation_system_charges
    
    def descriptor_atoms(self,molecules):
        """
        Yields descriptors for the input molecules based on the choosen descriptor.
        Descriptor uses the hyperparameter dictionary.
        
        
        Parameters
        ----------
        molecules : list
            list of ase atoms objects
            
        Returns
        -------
        descriptor : list 
            elements of the list depend on the choosen descriptor
            
        """ 
        descriptor = self.delta_descriptor(molecules)
        
        return descriptor

    def delta_descriptor(self,molecules):
        """
        Uses the atomic symbol as descriptor.
        
        Parameters
        ----------
        molecules : list
            list of ase atoms objects
        Returns
        -------
        descriptors : list
            descriptor is based on the atomic symbol
        """
        descriptors = [atom.symbol for mol in molecules for atom in mol]
        return descriptors


    def delta_kernel(self,mol_set1,mol_set2=None,kerneltype='general'):
        """
        Compares if two elements are equal and sets 1 if so, and 0 otherwise.
        
        Parameters
        ----------
        molsA : list 
            list of ase atoms objects
        compare_set : list
            list of ase atoms objects
            
        Returns
        -------
        K : matrix
            Kernel Matrix
        
        """
        if kerneltype == "training":
            desc1 = self.descriptor_atoms(self.training_set)
            dim1  = len(desc1)
            K = np.zeros((dim1,dim1))
            for iA in range(dim1):
                K[iA,iA] = 1.0
                for iB in range(iA):
                    if desc1[iA] == desc1[iB]:
                        K[iA,iB] = 1.0
                        K[iB,iA] = 1.0
        elif kerneltype == "predicting":
            desc2 = self.descriptor_atoms(self.training_set)
            desc1 = self.descriptor_atoms(mol_set1)
            dim1  = len(desc1)
            dim2  = len(desc2)
            K = np.zeros((dim1,dim2))
            
            for iA in range(dim1):
                for iB in range(dim2):
                    if desc1[iA] == desc2[iB]:
                        K[iA,iB] = 1.0
        elif kerneltype == "validation":
            desc2 = self.descriptor_atoms(self.training_set)
            desc1 = self.descriptor_atoms(self.validation_set)
            dim1  = len(desc1)
            dim2  = len(desc2)
            K = np.zeros((dim1,dim2))
            
            for iA in range(dim1):
                for iB in range(dim2):
                    if desc1[iA] == desc2[iB]:
                        K[iA,iB] = 1.0
        elif kerneltype == "general":
            if mol_set2 == None:
                desc1 = self.descriptor_atoms(mol_set1)
                dim1  = len(desc1)
                K = np.zeros((dim1,dim1))
                for iA in range(dim1):
                    K[iA,iA] = 1.0
                    for iB in range(iA):
                        if desc1[iA] == desc1[iB]:
                            K[iA,iB] = 1.0
                            K[iB,iA] = 1.0
            else:
                desc1 = self.descriptor_atoms(mol_set1)
                desc2 = self.descriptor_atoms(mol_set2)
                dim1  = len(desc1)
                dim2  = len(desc2)
                K = np.zeros((dim1,dim2))
                
                for iA in range(dim1):
                    for iB in range(dim2):
                        if desc1[iA] == desc2[iB]:
                            K[iA,iB] = 1.0
        return K


    def kernel_matrix(self,mol_set1=None,mol_set2=None,kerneltype='general'):
        """
        Returns the Kernel matrix for the specified kernel.
        
        Parameters
        ----------
        training_set : list
            list of ase atoms objects
        predict_set: list
            list of ase atoms objects
        
        Returns
        -------
        K : matrix
            Kernel Matrix
        """
        K = self.delta_kernel(mol_set1=mol_set1,mol_set2 = mol_set2,kerneltype=kerneltype)
        return K


    def _num_gradient(self,mol,desc_train,h=0.001,direction=0,iatom=0):
        tmpmol = mol.copy()
        pos = tmpmol.get_positions()/Bohr
        #pos /= Bohr
        pos[iatom][direction] += h
        tmpmol.set_positions(pos*Bohr)
        Kplus = self._soap_kernel_descriptors([tmpmol],desc_train)
        pos[iatom][direction] += -2.0*h
        tmpmol.set_positions(pos*Bohr)
        Kminus = self._soap_kernel_descriptors([tmpmol],desc_train)
        pos[iatom][direction] += h
        tmpmol.set_positions(pos*Bohr)
        return (Kplus-Kminus)/(2.0*h)
    
    def _soap_kernel_descriptors(self,ref_mol,descriptors_training,zeta=2):
        """
        The SOAP kernel between two atoms object and precomputed descriptors (for gradient code)
        
        Parameters
        ----------
        molA : obj
            ase atoms objects
        zeta : float
            exponent
        descriptors_training : list
            precomputed SOAP descriptors
            
        Returns
        -------
        K : matrix
            Kernel Matrix. Matrix elements range from 0 to 1
        """
        if self.multi_SOAP:
            deltas = self.descriptor_dict['deltas']
            desc1 = self.descriptor_atoms(ref_mol)
            dim1  = len(desc1)
            dim2  = len(descriptors_training)
            K_multi = []
            for soap_ind,soap_des in enumerate(desc1): 
                k = np.matmul(desc1[soap_ind],np.transpose(descriptors_training[soap_ind]))**zeta
                K_multi.append(np.multiply(deltas[soap_ind],k))
            K = 0
            for i in K_multi:
                K=np.add(i,K)
            return K
        else:
            desc1 = self.descriptor_atoms(ref_mol)
            #desc2 = self.descriptor_atoms(compare_set)
            #dim1  = len(desc1)
            dim2  = len(descriptors_training)
            K = np.matmul(desc1,np.transpose(descriptors_training))**zeta
            return K

    def _an_kernel_gradient(self,ref_mol,zeta=2):
        dKdr =  self._num_kernel_gradient(ref_mol,self.training_set)
        return dKdr