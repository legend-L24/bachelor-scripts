import numpy as np
from ase.io import *
from ase.data import covalent_radii
from ase.units import Bohr,Hartree
from kqeq.qeq import charge_eq
from kqeq.funct import _block_diag_rect,_block_diag,_get_R, get_energies, get_charges

from sklearn.metrics import mean_squared_error

class kernel_qeq():
    """Class for Kernel Charge Equilibration Models

    Parameters
    ----------
        Kernel: obj
            Kernel object with the function kernel_matrixQ
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
                include_LJ=False, radius_type='qeq', sparse=False, sparse_count=None, 
                sparse_method="CUR",hard_lib=None):
        
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
        self.K_train     = self.Kernel.kernel_matrix(kerneltype='training') #Kernelmatrix(self,self.training_set)
        self.A_train, self.X_train, self.O_train = self._build_A_X_O(mols=self.Kernel.training_set,systems_charges=self.Kernel.training_system_charges)
        self.R = None
        self.mu_ref = None
        
        self.q_ref = None
        if self.Kernel.validation_set:
            self.A_val, self.X_val, self.O_val = self._build_A_X_O(self.Kernel.validation_set,self.Kernel.validation_system_charges) 
            self.R_val = None
            self.mu_ref_val = None
            self.q_ref_val = None
            

        # Sparsification part, implemented FPS and CUR selection for representative set
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




    def update_kernel_matrix(self,New_Kernel):
        """Updates kernel matrix with a provided new kernel object.
 
        Parameters
        ----------
            New_Kernel : obj
                Kernel object with function kernel_matrix

        """
        self.Kernel     = New_Kernel
        self.K_train     = self.Kernel.kernel_matrix(kerneltype='training') #Kernelmatrix(self,self.training_set)
        self.A_train, self.X_train, self.O_train = self._build_A_X_O(mols=self.Kernel.training_set,systems_charges=self.Kernel.training_system_charges)
        self.R = None
        self.mu_ref = None
        self.q_ref = None
        if self.Kernel.validation_set:
            self.A_val, self.X_val, self.O_val = self._build_A_X_O(self.Kernel.validation_set,self.Kernel.validation_system_charges) 
            self.R_val = None
            self.mu_ref_val = None
            self.q_ref_val = None
            

        # Sparsification part, implemented FPS and CUR selection for representative set
        if self.sparse == True:
            if self.sparse_method == "FPS":
                self.rep_pos = self.farthest_point(self.K_train,self.sparse_count,1)
                self.Kernel.create_representationFPS(self.rep_pos)
            elif self.sparse_method == "CUR":
                self.sparse_points = self.Kernel.create_representationCUR(self.sparse_count)
            else:
                print("Specify sparse method!")
            self.K_nm, self.K_mm = self.Kernel.kernel_matrix(kerneltype='training')
##############################################################################################################################
################################### THIS IS CODE FOR KQEQ AND GAP ############################################################                                
##############################################################################################################################
    def train_GAPkQEq(self, lambda_reg = 0.1, lambda_charges = 0.001, iter_charges = 3, 
                    lambda_charges_min = 1e-6, charge_keyword='initial_charges', energy_keyword = "energy",
                    atom_energy = None):
        '''
        A is matrix from kQEq, so inverse hardness (A_up is zero matrix from GAP to this matrix)
        A_ones is matrix with ones and zeros which goes together with Q matrix from kQEq for energy
        At this moment only on additional target is allowed
        '''
        oldRMSE = 10000000
        K_nm = self.K_nm
        K_mm = self.K_mm        
        fullK_nm = np.kron(np.eye(2),K_nm)
        fullK_mm = np.kron(np.eye(2),K_mm)
        E_ref = get_energies(mols=self.Kernel.training_set,atom_energy = atom_energy, energy_keyword = energy_keyword)
        A, X, O, R, Xback = self._build_prediction(self.Kernel.training_set, self.Kernel.training_system_charges)
        if iter_charges > 0:
            T, ref = self._process_data("charges",self.Kernel.training_set,charge_keyword)
            nMol = len(self.Kernel.training_set)
            A_up = np.zeros((A.shape[0],A.shape[1]-nMol))
            fullA = np.hstack((A_up,A))
            X_up = np.zeros((X.shape[0]-nMol,X.shape[1]))
            fullX_up = np.hstack((X_up,X_up))
            fullX_temo = np.hstack((np.zeros((X.shape[0],X.shape[1])),X))
            fullX = np.vstack((fullX_up,fullX_temo))
            O_up = np.zeros(K_nm.shape[0])
            fullO = np.concatenate([O_up,O])
            target = None
            charges_part1 = np.linalg.multi_dot([fullK_nm.T,fullX.T,fullA.T,T.T])
            charges_part2 = np.linalg.multi_dot([T,fullA,fullO])
            charges_part3 = np.linalg.multi_dot([T,fullA,fullX,fullK_nm])
            #LambdaTarget = np.linalg.inv(lambda_target*np.eye(len(ref)))
        

        LambdaE = np.linalg.inv(lambda_reg*np.eye(len(E_ref)))
        A_ones = np.zeros((len(self.Kernel.training_set_descriptors),len(self.Kernel.training_set)))
        q_0 = get_charges(self.Kernel.training_set,charge_keyword=charge_keyword)
        Cs, Qmatrix = self._process_ref_energy(q_0)
        count_c = 0
        count_r = 0
        for mol in self.Kernel.training_set:
            for at in mol:
                A_ones[count_r,count_c] = 1
                count_r += 1
            count_c += 1
        train_en = True
        train_charges = False
        count_iter = 0
        bad = 0
        while train_en:
            M = np.vstack((A_ones,Qmatrix))
            up = np.linalg.multi_dot([fullK_nm.T,M,LambdaE,Cs]) - np.linalg.multi_dot([fullK_nm.T,M,LambdaE,E_ref])
            down = -np.linalg.multi_dot([fullK_nm.T,M, LambdaE,M.T,fullK_nm])
            if train_charges:
            
                Lambda_charges = ((1/lambda_charges)*np.eye(len(ref)))
                print("Target regularization:", lambda_charges)
                up += np.linalg.multi_dot([charges_part1,Lambda_charges,ref])
                up += np.linalg.multi_dot([charges_part1,Lambda_charges,charges_part2])
                down -= np.linalg.multi_dot([charges_part1,Lambda_charges,charges_part3])
                lambda_charges = lambda_charges/2
                if lambda_charges < lambda_charges_min:
                    lambda_charges = lambda_charges_min
                
            down -= fullK_mm

            weights = np.linalg.solve(down,up)
            WkQEQ = weights[K_nm.shape[1]:] 
            eneg  = np.matmul(K_nm,WkQEQ)
            eneg_tot = (np.matmul(X,np.transpose(eneg)))
            eneg_tot = eneg_tot + O
            charge_temp = np.matmul(A,-eneg_tot)
            q_new = np.matmul(Xback,charge_temp)
            ref = q_new
            
            #print(lambda_target)

            Cs_temp, Qmatrix_temp = self._process_ref_energy(q_new)
            M_temp = np.vstack((A_ones,Qmatrix))
            whole_energy =  np.linalg.multi_dot([M_temp.T,fullK_nm,weights])
            
            whole_energy = Cs_temp + whole_energy
            MSE = np.square(np.subtract(E_ref,whole_energy)).mean()  
            RMSE_E = np.sqrt(MSE)
            print("New RMSE:",RMSE_E)
            print("Old RMSE:",oldRMSE)
            print("Bad RMSEs in row:", bad)
            if RMSE_E > oldRMSE and train_charges == False:
                train_charges = True
                print("Charges added")
            elif RMSE_E > oldRMSE and train_charges == True and bad >= iter_charges:
                train_en = False
            elif RMSE_E > oldRMSE and train_charges == True and bad < iter_charges:
                bad += 1
                Cs = Cs_temp
                Qmatrix = Qmatrix_temp
                ref = q_new
            else:
                bad = 0
                final_weight = weights
                ref = q_new
                oldRMSE = RMSE_E
                Cs = Cs_temp
                Qmatrix = Qmatrix_temp
            print("Sample of charges:",q_new[-10:])
            print(f"iteration {count_iter} is done")
            count_iter += 1
        self.weights = final_weight 



    def calculate_GAPkQEq(self, mol, charge=0):
        """Calculates charges, energy, forces and dipole_vector for single atoms object. This is usable with trainEnergy function. After training electronegativities do not correspond with charges obtained via qe.calc_Charges()

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
        K = self.Kernel.kernel_matrix(mol_set1 = [mol],kerneltype="predicting")
        gapW = self.weights[:K.shape[1]]
        kQEqW = self.weights[K.shape[1]:]
        energyGAP = np.matmul(K,gapW)
        #print("GAP: ", np.sum(energyGAP)*Hartree)
        eneg   = np.matmul(K,kQEqW)
        
        qe     = charge_eq(mol,Qtot=charge,eneg=eneg,scale_atsize=self.scale_atsize,DoLJ=self.LJ,radius_type=self.radius_type, hard=self.hard_lib)
        qe.calc_Charges()
        qe.calc_Eqeq()
        mu = np.matmul(_get_R(mol),qe.q)
        energy  = qe.E 
        #print("kQEq:", energy*Hartree)
        #print("GAP:",np.sum(energyGAP)*Hartree)
        #print("kQEq:", energy*Hartree)
        return np.sum(energyGAP)*Hartree + energy*Hartree, qe.q, mu#, np.sum(energyGAP)*Hartree, energy*Hartree


    def trainEnCharge(self, lambda_reg = 0.1, lambda_target = 0.001, iter_charges = 5,
                   lambda_target_min = 1e-6, charge_keyword='initial_charges', 
                    energy_keyword = "energy",atom_energy = None):
        """
        this is a new function, which is the same as the latest GAP + kQEq training, 
        where only target is charge and energy with ingrising weight putted on charges. Probably delete later. 
        """
        T, ref = self._process_data("charges",self.Kernel.training_set,charge_keyword)
        oldRMSE = 10000
        E_ref = get_energies(mols=self.Kernel.training_set,atom_energy = atom_energy, energy_keyword = energy_keyword)
        A = self._build_A(self.Kernel.training_set)
        K_nm = self.K_nm
        K_mm = self.K_mm
        q_0 = get_charges(self.Kernel.training_set,charge_keyword=charge_keyword)
        q_0 = q_0.astype(float)
        Cs, Qmatrix = self._process_ref_energy(q_0)
        LambdaE = np.linalg.inv(lambda_reg*np.eye(len(E_ref)))
        #T, ref = self._process_data("charges",self.Kernel.training_set,charge_keyword)
        train_en = True
        train_charges = False
        count_iter = 0
        bad = 0
        while train_en:
            up = np.linalg.multi_dot([K_nm.T, Qmatrix, LambdaE,Cs])
            up = up - np.linalg.multi_dot([K_nm.T, Qmatrix, LambdaE,E_ref])
            down = np.linalg.multi_dot([K_nm.T,Qmatrix,LambdaE,Qmatrix.T,K_nm])
            
            if train_charges:
                up_temp, down_temp = self._train_up_down(T,ref,lambda_target)
                up += up_temp
                down += down_temp

            down += K_mm  
            weights = -np.linalg.solve(down,up)
            eneg_new  = np.matmul(K_nm,weights)
            q_new = np.matmul(A,-eneg_new)

            
            lambda_target = lambda_target/2
            if lambda_target < lambda_target_min:
                lambda_target = lambda_target_min
            Cs_temp, Qmatrix_temp = self._process_ref_energy(q_new)
            whole_energy =  np.linalg.multi_dot([Qmatrix_temp.T,K_nm,weights])
            whole_energy = Cs_temp + whole_energy
            MSE = np.square(np.subtract(E_ref,whole_energy)).mean()  
            RMSE_E = np.sqrt(MSE)
            # print("New RMSE:",RMSE_E)
            # print("Old RMSE:",oldRMSE)
            # print("Bad RMSEs in row:", bad)
            if count_iter == 50:
                train_en = False
            elif RMSE_E > oldRMSE and train_charges == False:
                train_charges = True
                # print("Charges added")
            elif RMSE_E > oldRMSE and train_charges == True and bad >= iter_charges:
                train_en = False
            elif RMSE_E > oldRMSE and train_charges == True and bad < iter_charges:
                bad += 1
                Cs = Cs_temp
                Qmatrix = Qmatrix_temp
                ref = q_new
            else:
                bad = 0
                final_weight = weights
                print("the final weight has been defined")
                ref = q_new
                oldRMSE = RMSE_E
                Cs = Cs_temp
                Qmatrix = Qmatrix_temp
            # print("Sample of charges:",q_new[-10:])
            # print(f"iteration {count_iter} is done")
            count_iter += 1
        self.weights = final_weight

    def trainEn(self, lambda_reg =0.1, Niter = 5, targets = [], 
                lambda_targets = [], charge_keyword='initial_charges',
                energy_keyword = "energy",atom_energy = None):
        """
        iterative training for fitting energy with regularization for electronegativities. Same as  trainEnQ2, but charges are deleted
        lambdaq and lambdaE has to be changed in the code itself, Niter is number of iterations
        Kw is weight for elektronegativties

        Need repair for charge target
        """
        if type(targets) == str:
            targets = [targets]
            lambda_targets = [lambda_targets]
        E_ref = get_energies(mols=self.Kernel.training_set,atom_energy = atom_energy, energy_keyword = energy_keyword)
        A = self._build_A(self.Kernel.training_set)
        K_nm = self.K_nm
        K_mm = self.K_mm
        q_0 = get_charges(self.Kernel.training_set,charge_keyword=charge_keyword)
        Cs, Qmatrix = self._process_ref_energy(q_0)
        #LambdaK = np.linalg.inv(lambda_reg*np.eye(K_nm.shape[0]))
        LambdaE = np.linalg.inv(lambda_reg*np.eye(len(E_ref)))
        for a in range(Niter):
            up = np.linalg.multi_dot([K_nm.T, Qmatrix, LambdaE,Cs])
            up = up - np.linalg.multi_dot([K_nm.T, Qmatrix, LambdaE,E_ref])
            down = np.linalg.multi_dot([K_nm.T,Qmatrix,LambdaE,Qmatrix.T,K_nm])
            for id,target in enumerate(targets):
                T, ref = self._process_data(target,self.Kernel.training_set,charge_keyword)
                up_temp, down_temp = self._train_up_down(T,ref,lambda_targets[id])
                up += up_temp
                down += down_temp
            #down +=  Kw*np.linalg.multi_dot([K_nm.T,LambdaK,K_nm])
            down += K_mm  
            weights = -np.linalg.solve(down,up)
            eneg_new  = np.matmul(K_nm,weights)
            q_new = np.matmul(A,-eneg_new)
            Cs, Qmatrix = self._process_ref_energy(q_new)
            print(q_new[-10:])
            print(f"iteration {a} is done")
        self.weights = weights

    def _process_ref_energy(self, charges):
        Qmatrix = np.zeros((len(self.Kernel.training_set_descriptors),len(self.Kernel.training_set)))
        Cs = []
        count_c = 0
        count_r = 0
        for mol in self.Kernel.training_set:
            q_temp = [] 
            for q in mol: 
                Qmatrix[count_r,count_c] = charges[count_r]
                q_temp.append(charges[count_r])
                count_r += 1
            count_c += 1
            #qe = charge_eq(mol,scale_atsize=self.scale_atsize,DoLJ=self.LJ,radius_type=self.radius_type, hard=self.hard_lib)
            qe = charge_eq(mol,Qtot = 0,scale_atsize=self.scale_atsize,DoLJ=self.LJ,radius_type=self.radius_type, hard=self.hard_lib)
            C = qe.compConst(q_temp)
            Cs.append(C)
        Cs = np.array(Cs)
        return Cs, Qmatrix

    def train(self, targets = "dipole", target_lambdas = 1, charge_keyword='initial_charges'):
        """Trains Kernel QEq model
 
        Parameters
        ----------
            lambda_reg : float
                Regularization hyperparameter
            target : string or list
                Current targets: "dipole", "charges" (can be combined e.g. ["dipole","charges"])
            target_weights : int or list
                Weights put on different targets. Use only when multiple targets are selected e.g. [0.8, 0.2]
            charge_keyword : string
                key for ASE .arrays dictionary for which charges should be used for "charges" training 
        """

        if type(targets) == str:
            targets = [targets]
            target_lambdas = [target_lambdas]
    
        up = np.zeros(self.K_mm.shape[0])
        down = np.zeros((self.K_mm.shape[0],self.K_mm.shape[0]))
        for id,target in enumerate(targets):
            T, ref = self._process_data(target,self.Kernel.training_set,charge_keyword)
            up_temp, down_temp = self._train_up_down(T,ref,target_lambdas[id])
            up += up_temp
            down += down_temp
        down = down+self.K_mm
        self.weights = -np.linalg.solve(down,up)
    

    def _train_up_down(self,T,ref,lambda_reg):
        A, K_nm, X, O = self.A_train, self.K_nm, self.X_train, self.O_train       
        Lambda = np.linalg.inv(lambda_reg*np.eye(len(ref)))
        up1 = np.linalg.multi_dot([K_nm.T,X.T,A.T,T.T,Lambda,ref])
        up2 = np.linalg.multi_dot([K_nm.T,X.T,A.T,T.T,Lambda,T,A,O])
        up = up1+up2
        down = np.linalg.multi_dot([K_nm.T,X.T,A.T,T.T,Lambda,T,A,X,K_nm])# + K_mm
        return up, down

        
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

    def _process_data(self,target, mols,charge_keyword):
        if target == "charges":
            ref_qs = []
            nAtoms = 0
            for m in mols:
                nAtoms += len(m)
            count_rowBack = 0
            count_colBack = 0
            Xback = np.zeros((nAtoms,nAtoms+len(mols)))
            for mol in mols:
                ref_qs.extend(mol.arrays[charge_keyword].astype(float))
                for at in mol:
                    Xback[count_rowBack,count_colBack] = 1
                    count_colBack += 1
                    count_rowBack += 1
                count_colBack += 1
                refdata = np.array(ref_qs)
            return Xback, refdata
        elif target == "dipole":
            all_rs = []
            ref_mus = []
            dimR = 0
            for mol in mols:
                ref_mus.extend(mol.info['dipole_vector'])
                all_rs.append(_get_R(mol))
                dimR += len(mol)+1

            R = _block_diag_rect(all_rs,dimR,len(all_rs)*3)
            refdata = np.array(ref_mus)
            return R, refdata



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
    

    def _farthest_point(self,K_matrix,n_struc,seeds):
        """Farthest point sampling algorithm for sparsification

        Parameters
        ----------
            K_matrix : array 
                Kernel matrix
            n_struc : int
                Number of structures for representative set
        
        Returns
        -------
            samples : array
                Representative set
        """   

        dia = K_matrix.diagonal()
        dist_matrix_squared = dia[:, None] + dia[None, :] - 2*K_matrix
        control_value = -1e-8
        control = np.nonzero(dist_matrix_squared < control_value)
        assert len(control[0]) == 0, f'Squared distance matrix element < {control_value}\n  \
                    Indices are : {np.argwhere(dist_matrix_squared<control_value)}.\n  \
                    Please check your code something went wrong!'

        dist_matrix_squared[control] = 0
        dist_matrix = np.sqrt(dist_matrix_squared)
        if isinstance(seeds, int):
            israise = True if n_struc <= 1 else False

            samples = [seeds]  # storage for furthest samples
            samples.append(np.argmax(dist_matrix[seeds]))  # find furthest from input sample

        elif isinstance(seeds, (list, np.ndarray)):
            israise = True if n_struc <= len(seeds) else False

            samples = [sample for sample in seeds]

        if israise:
            raise ValueError('`number` can not be smaller than specied in `seeds`')


        for idx in range(n_struc-len(samples)):
            samples_rem = np.delete(np.arange(len(dist_matrix)), samples)  # get indices of no selected samples

            dists = dist_matrix[samples][:, samples_rem]  # slice distances for selected samples to remaining samples

            dists_min = np.min(dists, axis=0)  # for each remaining sample find closest distance to already selected sample
            sample_furthest = np.argmax(dists_min)  # select the remaining sample furthest to all selected samples

            samples.append(samples_rem[sample_furthest])
        return samples


    def save_kQEq(self):
        np.save("weights.npy", self.weights)
        np.save("sparseSOAP.npy", self.Kernel.representing_set_descriptors)

        #np.savetxt("weights.kqeq", self.weights, delimiter="\n")
        #np.savetxt("sparseSOAP.kqeq", self.Kernel.representing_set_descriptors, delimiter="\n")
    
    def load_kQEq(self):
        self.weights = np.load("weights.npy")
        self.Kernel.representing_set_descriptors = np.load("sparseSOAP.npy")