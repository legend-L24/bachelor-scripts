
import numpy as np
from ase.units import Hartree,Bohr
from ase.data import covalent_radii, chemical_symbols
from scipy.special import erf, erfc
from kqeq.data import uff_xi_vdw, uff_Di_vdw, uff_Xi_qeq, uff_radius_qeq, initial_charge
from scipy.spatial.distance import cdist, pdist


#               H             Li            C             N             O
#eneg_default = {1:4.5/Hartree,3:4.5/Hartree,6:5.3/Hartree,7:6.9/Hartree,8:8.7/Hartree}

kcal2hartree = 0.0015936010974213599

class charge_eq():
    def __init__(self,mol,Qtot=0.0,q=None,eneg=[],hard=None,atsize=None,scale_atsize=1.0,DoLJ=False,radius_type='qeq'):
        """
        initialization and precomputation
        """
        # inputs and defaults
        self.xyz = mol.get_positions()/Bohr
        self.atoms = mol.get_atomic_numbers()
        self.nAt = len(self.atoms) 
        self.Qtot = Qtot
        self.DoLJ = DoLJ

        #this part added by ytli to do in peridic system
        self.pbc = mol.pbc.any()
        self.initial_charge = np.array([initial_charge[i] for i in mol.get_chemical_symbols()])
        if np.all(mol.cell.cellpar()[-3:]): 
            self.triclinic = 0
        else: 
            self.triclinic = 1
        self.accuracy = 0.000143996
        self.cutcol = 10
        self.qqrd2e = 14.3996
        self.cell = mol.get_cell()
        self.volume = mol.get_cell().volume

        #if q == None:
        if type(q) is np.ndarray:
            self.q = q
        elif q==None:
            self.q = np.zeros(self.nAt)
        else:
            raise ValueError("charge input is wrong in format(either None or array is allowed)")

        if len(eneg) == 0:
            eneg = []
            for El in self.atoms:
                eneg.append(uff_Xi_qeq[chemical_symbols[El]]/Hartree)
            self.eneg = np.array(eneg)
        else:
            self.eneg = eneg

        if hard == None:
            self.hard = np.zeros(self.nAt)
        else:
            hard_temp = []
            for El in self.atoms:
                hard_temp.append(hard[chemical_symbols[El]])
            self.hard = np.array(hard_temp)


        if atsize == None:
            atsize = []
            for El in self.atoms:
                if radius_type=='rcov':
                    atsize.append(covalent_radii[El]/Bohr) 
                elif radius_type=='qeq':
                    atsize.append(uff_radius_qeq[chemical_symbols[El]]/Bohr) 
            self.atsize = scale_atsize*np.array(atsize)
        else:
            self.atsize = atsize

        if self.DoLJ:
            xi_LJ = [] 
            Di_LJ = []
            for El in self.atoms:
                xi_LJ.append(uff_xi_vdw[chemical_symbols[El]]/Bohr)
                Di_LJ.append(uff_Di_vdw[chemical_symbols[El]]*kcal2hartree)
            self.xi_LJ = np.array(xi_LJ)
            self.Di_LJ = np.array(Di_LJ)
        
        # calc distance matrix
        self.calc_Rij()

        self.E = None
    # a loss function for ewald summation
    def rms(self, km, prd, natom, q2, g_ewald):
        return (2*q2*g_ewald/prd)*np.sqrt(1/np.pi/km/natom)*np.exp(-np.pi*np.pi*km*km/((g_ewald*prd)**2))
    
    def calc_Rij(self):
        self.Rij = np.zeros((self.nAt,self.nAt))
        for i in range(self.nAt):
            for j in range(self.nAt):
                self.Rij[i,j] = self.calc_Rvec(i,j)

    def calc_Rvec(self,i,j):
        vi = self.xyz[i]
        vj = self.xyz[j]
        rij = np.subtract(vj,vi)
        r = np.linalg.norm(rij)
        return r


    def calc_dVdr(self,ga_ij,Rij):
        dVdr = 0.0
        if Rij>0.0:
            dVdr += 2.0*ga_ij*np.exp(-ga_ij**2 * Rij**2)/(np.sqrt(np.pi)*Rij)
            dVdr += -erf(ga_ij*Rij)/(Rij**2)
        return dVdr

    def calc_ELJ(self):
        etmp = 0.0
        for i in range(self.nAt):
            for j in range(i):
                rvec = self.xyz[i] - self.xyz[j]
                rij = np.linalg.norm(rvec)
                
                Dij = np.sqrt(self.Di_LJ[i]*self.Di_LJ[j])
                xij = np.sqrt(self.xi_LJ[i]*self.xi_LJ[j])
                xr6  = (xij/rij)**6
                xr12 = xr6**2
                etmp += Dij*(-2.*xr6+xr12)
                
                #Forces
                dVdr = 12.*(Dij/rij)*(xr6-xr12)
                self.f[i,:] += -dVdr/rij * rvec
                self.f[j,:] +=  dVdr/rij * rvec

        self.E += etmp
    
    def calc_gamma(self,a1,a2):
        return 1.0/np.sqrt(a1**2+a2**2)

    def calc_Vscreen(self,ga_ij,Rij):
        return erf(ga_ij*Rij)/Rij


    def calc_Eqeq_old(self):
        self.E = 0
        self.f = np.zeros((self.nAt,3))
        # Hardness and Electronegativity contribution
        E = np.sum(self.q*self.eneg) + 0.5*np.sum(self.q**2 *self.hard)

        # Set up matrices for Coulomb
        Rij = np.sqrt(np.sum((self.xyz[:,np.newaxis,:]-self.xyz[np.newaxis,:,:])**2,axis=-1)) # see scipy distance matrix
        eye = np.eye(Rij.shape[0])
        Rij = Rij + eye
        antieye = np.ones_like(Rij) - eye
        ga_ij = 1.0/np.sqrt(self.atsize[:,np.newaxis]**2 + self.atsize[np.newaxis,:]**2)
        Vscreen = erf(ga_ij*Rij)/Rij
        Vscreen = Vscreen * antieye
        E = E + 0.5*self.q@Vscreen@self.q

        # Local Coulomb self-energy
        ga_ii = np.diag(ga_ij)
        E = E + np.sum(0.5*(2.0*ga_ii/np.sqrt(np.pi))*self.q**2)
        self.E = E
        dVdr = 2.0*ga_ij*np.exp(-ga_ij**2 * Rij**2)/(np.sqrt(np.pi)*Rij)-erf(ga_ij*Rij)/(Rij**2)
        middle = self.q*dVdr/Rij
        for a in range((self.nAt)):
            self.f += self.q[a]*(self.xyz[a]-self.xyz)*middle[a][:, np.newaxis]
    
    def calc_Eqeq(self):
        H = self.get_H()
        E = np.sum(self.q*self.eneg)
        self.f = np.zeros((self.nAt,3))
        self.E = E + 0.5*self.q@H@self.q


    def calc_Eqeq_old(self):
        #print("calculating electrostatic energy")
        #print(self.q)
        self.E = 0.0
        self.f = np.zeros((self.nAt,3))
        etmp = 0.0
        for i in range(self.nAt):
            qi = self.q[i]
            ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])
        
            etmp += self.eneg[i]*qi
            etmp += 0.5*(self.hard[i]+2.0*ga_ii/np.sqrt(np.pi))*qi**2
            #print(self.hard[i]+2.0*ga_ii/np.sqrt(np.pi))
            #etmp += (2.0*ga_ii/np.sqrt(np.pi))*qi**2
            for j in range(i):
                qj = self.q[j] 
                ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                #print("ga_ij:", ga_ij)
                rij = self.Rij[i,j]
                etmp += qi*qj*self.calc_Vscreen(ga_ij,rij)
                #Forces
                rvec = self.xyz[i] - self.xyz[j]
                dVdr = self.calc_dVdr(ga_ij,rij)
                self.f[i,:] += -qi*qj*dVdr/rij * rvec
                self.f[j,:] +=  qi*qj*dVdr/rij * rvec
        self.E = etmp
    

    def calc_ChemPot(self,i):
        #A = self.atoms[i]
        #para = self.par[A]
        qi = self.q[i]
        Xi = self.eneg[i]
        for j in range(self.nAt):
            qj = self.q[j]
            ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
            rij = self.Rij[i,j]
            Xi += self.calc_Vscreen(i,j)*qj
        return Xi

    def get_Abar(self):
        if self.pbc:
            print("this is pbc examples")
            q2 = self.qqrd2e*np.sum(np.square(self.initial_charge))
            cell_vec = self.cell
            natom = self.nAt
            cutoff = self.cutcol
            xprd = cell_vec[0][0]
            yprd = cell_vec[1][1]
            zprd = cell_vec[2][2]
            g_ewald = self.accuracy*np.sqrt(natom*cutoff*xprd*yprd*zprd) / (2.0*q2)
            if (g_ewald >= 1.0): g_ewald = (1.35 - 0.15*np.log(self.accuracy))/cutoff
            else: g_ewald = np.sqrt(-np.log(g_ewald)) / cutoff
            la_vec = self.cell[:]
            #generate k-space:
            a = np.copy(la_vec)
            a[[0,1], :] = a[[1,0], :]
            a[[1,2], :] = a[[2,1], :]

            b = np.copy(la_vec)
            b[[1,2], :] = b[[2,1], :]
            b[[0,1], :] = b[[1,0], :]
            k_vec = 2*np.pi*np.cross(a,b)/self.volume

            pos_ls = []
            car_pos = self.xyz*Bohr
            for i in car_pos:
                for j in car_pos:
                    pos_ls.append(i-j)
            pos_vec = np.array(pos_ls)
            # set up the range in real and recipical space
            # check the k-space limit

            #accuracy = 0.000143996 # 1e-4
            
            
            kxmax = kymax = kzmax = 1
            err = self.rms(kxmax,xprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kxmax+=1
                err = self.rms(kxmax,xprd,natom,q2,g_ewald)
            err = self.rms(kymax,yprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kymax+=1
                err = self.rms(kymax,yprd,natom,q2,g_ewald)
            err = self.rms(kzmax,zprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kzmax+=1
                err = self.rms(kzmax,zprd,natom,q2,g_ewald)
                
            kmax_vec = np.array([kxmax,kymax,kzmax])
            kprd = np.array([1/xprd, 1/yprd, 1/zprd]) * 2*np.pi
            gsqmx = np.max(kmax_vec*kprd)

                
            if self.triclinic:
                h_vec = np.array([cell_vec[0,0], cell_vec[1,1], cell_vec[2,2], abs(cell_vec[2,1]), abs(cell_vec[2,0]), abs(cell_vec[1,0])])  #shape matrix in Voigt ordering; Voigt = xx,yy,zz,yz,xz,xy
                lamda = np.array([kxmax, kymax, kzmax])/np.array([xprd, yprd, zprd])
                v_tmp = np.array([0,0,0])
                v_tmp[0] = h_vec[0]*lamda[0]
                v_tmp[1] = h_vec[5]*lamda[0] + h_vec[1]*lamda[1]
                v_tmp[2] = h_vec[4]*lamda[0] + h_vec[3]*lamda[1] + h_vec[2]*lamda[2]
                v_tmp = v_tmp.astype(np.int32)
                kmax_vec = np.max(np.vstack((np.array([1,1,1]), v_tmp)), axis=0)

            kmax = np.max(kmax_vec)
            gsqmx *= 1.00001
            real_limit = [0,0,0]
            for i in range(0,3):
                dis_vec = la_vec[i]
                while 1: 
                    real_limit[i] += 1
                    dis_pos = real_limit[i]*dis_vec + car_pos
                    dist_mat = cdist(dis_pos,car_pos)
                    if np.min(dist_mat) > self.cutcol: break
            dis_ls = []
            kvec_ls = []

            #(k,0,0), (0,l,0), (0,0,m)
            for i in range(1, kmax+1):
                kvec_ls.append([i,0,0])
                kvec_ls.append([0,i,0])
                kvec_ls.append([0,0,i])
            # 1 = (k,l,0), 2 = (k,-l,0)
            for k in range(1, kmax_vec[0]+1):
                for l in range(1, kmax_vec[1]+1):
                    kvec_ls.append([k,l,0])
                    kvec_ls.append([k,-l,0])
            # 1 = (0,l,m), 2 = (0,l,-m)
            for m in range(1, kmax_vec[2]+1):
                for l in range(1, kmax_vec[1]+1):
                    kvec_ls.append([0,l,m])
                    kvec_ls.append([0,l,-m])
            # 1 = (k,0,m), 2 = (k,0,-m)
            for k in range(1, kmax_vec[0]+1):
                for m in range(1, kmax_vec[2]+1):
                    kvec_ls.append([k,0,m])
                    kvec_ls.append([k,0,-m])
            # 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
            for k in range(1, kmax_vec[0]+1):
                for l in range(1, kmax_vec[1]+1):
                    for m in range(1, kmax_vec[2]+1):
                        kvec_ls.append([k,l,m])
                        kvec_ls.append([k,l,-m])
                        kvec_ls.append([k,-l,m])
                        kvec_ls.append([k,-l,-m])

            kvec = np.array(kvec_ls)
            kvec = np.dot(kvec,k_vec)
            kvec = np.hstack((kvec,np.linalg.norm(kvec,axis=1).reshape(kvec.shape[0],1)))
            kvec = kvec[kvec[:,-1]<gsqmx]

            phase_vec = np.dot(kvec[:,:3], pos_vec.T)

            # k-space part
            sqk = np.square(kvec[:,-1])
            ug_vec = 8*np.pi*np.exp(-0.25*sqk/g_ewald/g_ewald)/self.volume/sqk # the ug_vec double since we merge (h,k,l) and (-h,-k,-l)
            sfactor = np.cos(phase_vec)
            A_bar_reci = np.dot(sfactor.T,ug_vec).reshape(natom, natom)

            # self energy
            A_bar_reci -= 2* np.eye(natom,natom)*g_ewald/np.sqrt(np.pi)

            for i in range(-real_limit[0], real_limit[0]+1):
                for j in range(-real_limit[1], real_limit[1]+1):
                    for k in range(-real_limit[2], real_limit[2]+1):
                        dis_vec = np.dot(np.array([i,j,k]),la_vec)
                        dis_pos = dis_vec + car_pos
                        dist_mat = cdist(dis_pos, car_pos)
                        if np.min(dist_mat) > cutoff: continue
                        pbc_mat = np.dot(np.array([i,j,k]),la_vec)
                        dis_ls.append(car_pos+pbc_mat)
                        #dis_ls.append(dis_pos)
            dis_num = len(dis_ls)
            diff_pos = np.concatenate(dis_ls)
            size_vec = np.transpose(np.matrix(self.atsize*Bohr))
            ga_ii = np.square(size_vec)
            gama_mat = np.sqrt(ga_ii+ga_ii.T)
            A_bar = np.zeros((natom,natom))
            # case 1: i != j

            dist_mat = cdist(car_pos, diff_pos)
            dist_mat[dist_mat>cutoff] = 0
            A_bar_all = np.dot(np.array([i,j,k]),la_vec)
            gama_mat = np.concatenate([gama_mat]*dis_num,axis=1)
            A_bar_all = np.divide(erfc(g_ewald*dist_mat)-erfc(dist_mat/np.sqrt(2)/gama_mat), dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            for ii in range(0, dis_num):
                A_bar += A_bar_all[:, ii*natom: (ii+1)*natom]

            #A_bar = np.divide(erfc(g_ewald*dist_mat), dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            #A_bar += np.divide(erfc(g_ewald*dist_mat)-erfc(dist_mat/np.sqrt(2)/gama_mat) , dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            # case 2: i == j
            A_bar = A_bar + np.eye(self.nAt, self.nAt)/np.sqrt(np.pi)/size_vec
            A_bar = np.array(A_bar)
            A = A_bar + A_bar_reci
            A = np.row_stack((A, np.ones(self.nAt)))
            A = np.column_stack((A, np.ones(self.nAt+1)))
            A[-1, -1] = 0
        else:
            A = np.zeros((self.nAt+1,self.nAt+1))
            # set up coefficient matrix
            for i in range(self.nAt):
                #para = self.par[ElA]
                ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])
                A[i,-1] = 1.0
                for j in range(self.nAt):
                    #parb = self.par[ElB]
                    if i==j:
                        A[i,j] = self.hard[i]+2.0*ga_ii/np.sqrt(np.pi)
                    else:
                        ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                        rij = self.Rij[i,j]
                        A[i,j] = self.calc_Vscreen(ga_ij,rij)
            for i in range(self.nAt):
                A[-1,i] = 1.0
            A[-1,-1] = 0.0
        Ainv = np.linalg.inv(A)
        return Ainv

    def get_H(self):
        A = np.zeros((self.nAt,self.nAt))
        if self.pbc:
            q2 = self.qqrd2e*np.sum(np.square(self.initial_charge))
            cell_vec = self.cell
            natom = self.nAt
            cutoff = self.cutcol
            xprd = cell_vec[0][0]
            yprd = cell_vec[1][1]
            zprd = cell_vec[2][2]
            g_ewald = self.accuracy*np.sqrt(natom*cutoff*xprd*yprd*zprd) / (2.0*q2)
            if (g_ewald >= 1.0): g_ewald = (1.35 - 0.15*np.log(self.accuracy))/cutoff
            else: g_ewald = np.sqrt(-np.log(g_ewald)) / cutoff
            la_vec = self.cell[:]
            #generate k-space:
            a = np.copy(la_vec)
            a[[0,1], :] = a[[1,0], :]
            a[[1,2], :] = a[[2,1], :]

            b = np.copy(la_vec)
            b[[1,2], :] = b[[2,1], :]
            b[[0,1], :] = b[[1,0], :]
            k_vec = 2*np.pi*np.cross(a,b)/self.volume

            pos_ls = []
            car_pos = self.xyz*Bohr
            for i in car_pos:
                for j in car_pos:
                    pos_ls.append(i-j)
            pos_vec = np.array(pos_ls)
            # set up the range in real and recipical space
            # check the k-space limit

            #accuracy = 0.000143996 # 1e-4
            

            kxmax = kymax = kzmax = 1
            err = self.rms(kxmax,xprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kxmax+=1
                err = self.rms(kxmax,xprd,natom,q2,g_ewald)
            err = self.rms(kymax,yprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kymax+=1
                err = self.rms(kymax,yprd,natom,q2,g_ewald)
            err = self.rms(kzmax,zprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kzmax+=1
                err = self.rms(kzmax,zprd,natom,q2,g_ewald)
                
            kmax_vec = np.array([kxmax,kymax,kzmax])
            kprd = np.array([1/xprd, 1/yprd, 1/zprd]) * 2*np.pi
            gsqmx = np.max(kmax_vec*kprd)

                
            if self.triclinic:
                h_vec = np.array([cell_vec[0,0], cell_vec[1,1], cell_vec[2,2], abs(cell_vec[2,1]), abs(cell_vec[2,0]), abs(cell_vec[1,0])])  #shape matrix in Voigt ordering; Voigt = xx,yy,zz,yz,xz,xy
                lamda = np.array([kxmax, kymax, kzmax])/np.array([xprd, yprd, zprd])
                v_tmp = np.array([0,0,0])
                v_tmp[0] = h_vec[0]*lamda[0]
                v_tmp[1] = h_vec[5]*lamda[0] + h_vec[1]*lamda[1]
                v_tmp[2] = h_vec[4]*lamda[0] + h_vec[3]*lamda[1] + h_vec[2]*lamda[2]
                v_tmp = v_tmp.astype(np.int32)
                kmax_vec = np.max(np.vstack((np.array([1,1,1]), v_tmp)), axis=0)

            kmax = np.max(kmax_vec)
            gsqmx *= 1.00001
            real_limit = [0,0,0]
            for i in range(0,3):
                dis_vec = la_vec[i]
                while 1: 
                    real_limit[i] += 1
                    dis_pos = real_limit[i]*dis_vec + car_pos
                    dist_mat = cdist(dis_pos,car_pos)
                    if np.min(dist_mat) > self.cutcol: break
            dis_ls = []
            kvec_ls = []

            #(k,0,0), (0,l,0), (0,0,m)
            for i in range(1, kmax+1):
                kvec_ls.append([i,0,0])
                kvec_ls.append([0,i,0])
                kvec_ls.append([0,0,i])
            # 1 = (k,l,0), 2 = (k,-l,0)
            for k in range(1, kmax_vec[0]+1):
                for l in range(1, kmax_vec[1]+1):
                    kvec_ls.append([k,l,0])
                    kvec_ls.append([k,-l,0])
            # 1 = (0,l,m), 2 = (0,l,-m)
            for m in range(1, kmax_vec[2]+1):
                for l in range(1, kmax_vec[1]+1):
                    kvec_ls.append([0,l,m])
                    kvec_ls.append([0,l,-m])
            # 1 = (k,0,m), 2 = (k,0,-m)
            for k in range(1, kmax_vec[0]+1):
                for m in range(1, kmax_vec[2]+1):
                    kvec_ls.append([k,0,m])
                    kvec_ls.append([k,0,-m])
            # 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
            for k in range(1, kmax_vec[0]+1):
                for l in range(1, kmax_vec[1]+1):
                    for m in range(1, kmax_vec[2]+1):
                        kvec_ls.append([k,l,m])
                        kvec_ls.append([k,l,-m])
                        kvec_ls.append([k,-l,m])
                        kvec_ls.append([k,-l,-m])

            kvec = np.array(kvec_ls)
            kvec = np.dot(kvec,k_vec)
            kvec = np.hstack((kvec,np.linalg.norm(kvec,axis=1).reshape(kvec.shape[0],1)))
            kvec = kvec[kvec[:,-1]<gsqmx]

            phase_vec = np.dot(kvec[:,:3], pos_vec.T)

            # k-space part
            sqk = np.square(kvec[:,-1])
            ug_vec = 8*np.pi*np.exp(-0.25*sqk/g_ewald/g_ewald)/self.volume/sqk # the ug_vec double since we merge (h,k,l) and (-h,-k,-l)
            sfactor = np.cos(phase_vec)
            A_bar_reci = np.dot(sfactor.T,ug_vec).reshape(natom, natom)

            # self energy
            A_bar_reci -= 2* np.eye(natom,natom)*g_ewald/np.sqrt(np.pi)

            for i in range(-real_limit[0], real_limit[0]+1):
                for j in range(-real_limit[1], real_limit[1]+1):
                    for k in range(-real_limit[2], real_limit[2]+1):
                        dis_vec = np.dot(np.array([i,j,k]),la_vec)
                        dis_pos = dis_vec + car_pos
                        dist_mat = cdist(dis_pos, car_pos)
                        if np.min(dist_mat) > cutoff: continue
                        pbc_mat = np.dot(np.array([i,j,k]),la_vec)
                        dis_ls.append(car_pos+pbc_mat)
                        #dis_ls.append(dis_pos)
            dis_num = len(dis_ls)
            diff_pos = np.concatenate(dis_ls)
            size_vec = np.transpose(np.matrix(self.atsize*Bohr))
            ga_ii = np.square(size_vec)
            gama_mat = np.sqrt(ga_ii+ga_ii.T)
            A_bar = np.zeros((natom,natom))
            # case 1: i != j

            dist_mat = cdist(car_pos, diff_pos)
            dist_mat[dist_mat>cutoff] = 0
            A_bar_all = np.dot(np.array([i,j,k]),la_vec)
            gama_mat = np.concatenate([gama_mat]*dis_num,axis=1)
            A_bar_all = np.divide(erfc(g_ewald*dist_mat)-erfc(dist_mat/np.sqrt(2)/gama_mat), dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            for ii in range(0, dis_num):
                A_bar += A_bar_all[:, ii*natom: (ii+1)*natom]

            #A_bar = np.divide(erfc(g_ewald*dist_mat), dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            #A_bar += np.divide(erfc(g_ewald*dist_mat)-erfc(dist_mat/np.sqrt(2)/gama_mat) , dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            # case 2: i == j
            A_bar = A_bar + np.eye(self.nAt, self.nAt)/np.sqrt(np.pi)/size_vec
            A_bar = np.array(A_bar)
            A = A_bar + A_bar_reci
            return A
        else: # set up coefficient matrix
            for i in range(self.nAt):
                #para = self.par[ElA]
                ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])
                A[i,-1] = 1.0
                for j in range(self.nAt):
                    #parb = self.par[ElB]
                    if i==j:
                        A[i,j] = self.hard[i]+2.0*ga_ii/np.sqrt(np.pi)
                        #print("this is {0} {1}: ".format(i,j), A[i,j])
                    else:
                        ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                        rij = self.Rij[i,j]
                        A[i,j] = self.calc_Vscreen(ga_ij,rij)
            return A

    def get_A(self):
        A = np.zeros((self.nAt+1,self.nAt+1))
        if self.pbc:
            print("this is pbc examples")
            q2 = self.qqrd2e*np.sum(np.square(self.initial_charge))
            cell_vec = self.cell
            natom = self.nAt
            cutoff = self.cutcol
            xprd = cell_vec[0][0]
            yprd = cell_vec[1][1]
            zprd = cell_vec[2][2]
            g_ewald = self.accuracy*np.sqrt(natom*cutoff*xprd*yprd*zprd) / (2.0*q2)
            if (g_ewald >= 1.0): g_ewald = (1.35 - 0.15*np.log(self.accuracy))/cutoff
            else: g_ewald = np.sqrt(-np.log(g_ewald)) / cutoff
            la_vec = self.cell[:]
            #generate k-space:
            a = np.copy(la_vec)
            a[[0,1], :] = a[[1,0], :]
            a[[1,2], :] = a[[2,1], :]

            b = np.copy(la_vec)
            b[[1,2], :] = b[[2,1], :]
            b[[0,1], :] = b[[1,0], :]
            k_vec = 2*np.pi*np.cross(a,b)/self.volume

            pos_ls = []
            car_pos = self.xyz*Bohr
            for i in car_pos:
                for j in car_pos:
                    pos_ls.append(i-j)
            pos_vec = np.array(pos_ls)
            # set up the range in real and recipical space
            # check the k-space limit

            #accuracy = 0.000143996 # 1e-4
            

            kxmax = kymax = kzmax = 1
            err = self.rms(kxmax,xprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kxmax+=1
                err = self.rms(kxmax,xprd,natom,q2,g_ewald)
            err = self.rms(kymax,yprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kymax+=1
                err = self.rms(kymax,yprd,natom,q2,g_ewald)
            err = self.rms(kzmax,zprd,natom,q2,g_ewald)
            while (err > self.accuracy):
                kzmax+=1
                err = self.rms(kzmax,zprd,natom,q2,g_ewald)
                
            kmax_vec = np.array([kxmax,kymax,kzmax])
            kprd = np.array([1/xprd, 1/yprd, 1/zprd]) * 2*np.pi
            gsqmx = np.max(kmax_vec*kprd)

                
            if self.triclinic:
                h_vec = np.array([cell_vec[0,0], cell_vec[1,1], cell_vec[2,2], abs(cell_vec[2,1]), abs(cell_vec[2,0]), abs(cell_vec[1,0])])  #shape matrix in Voigt ordering; Voigt = xx,yy,zz,yz,xz,xy
                lamda = np.array([kxmax, kymax, kzmax])/np.array([xprd, yprd, zprd])
                v_tmp = np.array([0,0,0])
                v_tmp[0] = h_vec[0]*lamda[0]
                v_tmp[1] = h_vec[5]*lamda[0] + h_vec[1]*lamda[1]
                v_tmp[2] = h_vec[4]*lamda[0] + h_vec[3]*lamda[1] + h_vec[2]*lamda[2]
                v_tmp = v_tmp.astype(np.int32)
                kmax_vec = np.max(np.vstack((np.array([1,1,1]), v_tmp)), axis=0)

            kmax = np.max(kmax_vec)
            gsqmx *= 1.00001
            real_limit = [0,0,0]
            for i in range(0,3):
                dis_vec = la_vec[i]
                while 1: 
                    real_limit[i] += 1
                    dis_pos = real_limit[i]*dis_vec + car_pos
                    dist_mat = cdist(dis_pos,car_pos)
                    if np.min(dist_mat) > self.cutcol: break
            dis_ls = []
            kvec_ls = []

            #(k,0,0), (0,l,0), (0,0,m)
            for i in range(1, kmax+1):
                kvec_ls.append([i,0,0])
                kvec_ls.append([0,i,0])
                kvec_ls.append([0,0,i])
            # 1 = (k,l,0), 2 = (k,-l,0)
            for k in range(1, kmax_vec[0]+1):
                for l in range(1, kmax_vec[1]+1):
                    kvec_ls.append([k,l,0])
                    kvec_ls.append([k,-l,0])
            # 1 = (0,l,m), 2 = (0,l,-m)
            for m in range(1, kmax_vec[2]+1):
                for l in range(1, kmax_vec[1]+1):
                    kvec_ls.append([0,l,m])
                    kvec_ls.append([0,l,-m])
            # 1 = (k,0,m), 2 = (k,0,-m)
            for k in range(1, kmax_vec[0]+1):
                for m in range(1, kmax_vec[2]+1):
                    kvec_ls.append([k,0,m])
                    kvec_ls.append([k,0,-m])
            # 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
            for k in range(1, kmax_vec[0]+1):
                for l in range(1, kmax_vec[1]+1):
                    for m in range(1, kmax_vec[2]+1):
                        kvec_ls.append([k,l,m])
                        kvec_ls.append([k,l,-m])
                        kvec_ls.append([k,-l,m])
                        kvec_ls.append([k,-l,-m])

            kvec = np.array(kvec_ls)
            kvec = np.dot(kvec,k_vec)
            kvec = np.hstack((kvec,np.linalg.norm(kvec,axis=1).reshape(kvec.shape[0],1)))
            kvec = kvec[kvec[:,-1]<gsqmx]

            phase_vec = np.dot(kvec[:,:3], pos_vec.T)

            # k-space part
            sqk = np.square(kvec[:,-1])
            ug_vec = 8*np.pi*np.exp(-0.25*sqk/g_ewald/g_ewald)/self.volume/sqk # the ug_vec double since we merge (h,k,l) and (-h,-k,-l)
            sfactor = np.cos(phase_vec)
            A_bar_reci = np.dot(sfactor.T,ug_vec).reshape(natom, natom)

            # self energy
            A_bar_reci -= 2* np.eye(natom,natom)*g_ewald/np.sqrt(np.pi)

            for i in range(-real_limit[0], real_limit[0]+1):
                for j in range(-real_limit[1], real_limit[1]+1):
                    for k in range(-real_limit[2], real_limit[2]+1):
                        dis_vec = np.dot(np.array([i,j,k]),la_vec)
                        dis_pos = dis_vec + car_pos
                        dist_mat = cdist(dis_pos, car_pos)
                        if np.min(dist_mat) > cutoff: continue
                        pbc_mat = np.dot(np.array([i,j,k]),la_vec)
                        dis_ls.append(car_pos+pbc_mat)
                        #dis_ls.append(dis_pos)
            dis_num = len(dis_ls)
            diff_pos = np.concatenate(dis_ls)
            size_vec = np.transpose(np.matrix(self.atsize*Bohr))
            ga_ii = np.square(size_vec)
            gama_mat = np.sqrt(ga_ii+ga_ii.T)
            A_bar = np.zeros((natom,natom))
            # case 1: i != j

            dist_mat = cdist(car_pos, diff_pos)
            dist_mat[dist_mat>cutoff] = 0
            A_bar_all = np.dot(np.array([i,j,k]),la_vec)
            gama_mat = np.concatenate([gama_mat]*dis_num,axis=1)
            A_bar_all = np.divide(erfc(g_ewald*dist_mat)-erfc(dist_mat/np.sqrt(2)/gama_mat), dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            for ii in range(0, dis_num):
                A_bar += A_bar_all[:, ii*natom: (ii+1)*natom]

            #A_bar = np.divide(erfc(g_ewald*dist_mat), dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            #A_bar += np.divide(erfc(g_ewald*dist_mat)-erfc(dist_mat/np.sqrt(2)/gama_mat) , dist_mat, out=np.zeros_like(dist_mat), where=dist_mat!=0)
            # case 2: i == j
            A_bar = A_bar + np.eye(self.nAt, self.nAt)/np.sqrt(np.pi)/size_vec
            A_bar = np.array(A_bar)
            A = A_bar + A_bar_reci
            A = np.row_stack((A, np.ones(self.nAt)))
            A = np.column_stack((A, np.ones(self.nAt+1)))
            A[-1, -1] = 0
            #self.A = A
            #return np.linalg.inv(A)
        else:
            # set up coefficient matrix
            for i in range(self.nAt):
                #para = self.par[ElA]
                ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])
                A[i,-1] = 1.0
                for j in range(self.nAt):
                    #parb = self.par[ElB]
                    if i==j:
                        A[i,j] = self.hard[i]+2.0*ga_ii/np.sqrt(np.pi)
                    else:
                        ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                        rij = self.Rij[i,j]
                        A[i,j] = self.calc_Vscreen(ga_ij,rij)
            for i in range(self.nAt):
                A[-1,i] = 1.0
            A[-1,-1] = 0.0
        Ainv = np.linalg.inv(A)
        return Ainv[:-1,:-1]

    def calc_Charges(self):
        A = np.zeros((self.nAt+1,self.nAt+1))
        #q = np.zeros((self.nAt+1))
        C = np.zeros(self.nAt+1)
        # set up solutions vector
        for i in range(self.nAt):
            #para = self.par[ElA]
            C[i] = -self.eneg[i]
        C[-1] = self.Qtot 
        # set up coefficient matrix
        '''
        for i in range(self.nAt):
            #para = self.par[ElA]
            ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])
            A[i,-1] = 1.0
            for j in range(self.nAt):
                #parb = self.par[ElB]
                if i==j:
                    A[i,j] = self.hard[i]+2.0*ga_ii/np.sqrt(np.pi)    
                else:
                    ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                    rij = self.Rij[i,j]
                    A[i,j] = self.calc_Vscreen(ga_ij,rij)
        for i in range(self.nAt):
            A[-1,i] = 1.0 
        A[-1,-1] = 0.0
        #print(A)
        Ainv = np.linalg.inv(A)
        '''
        Ainv = self.get_Abar()
        q = np.matmul(Ainv,C)
        #q = np.linalg.lstsq(A,C,rcond=None)
        self.q = q[:-1]
        #print(q)

    def compConst(self,charges):
        H = self.get_H()
        char_arr = np.array(charges)
        C = 0.5*char_arr@H@char_arr
        '''
        for i in range(self.nAt):
            qi = charges[i]
            ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])
            C += 0.5*(self.hard[i]+2.0*ga_ii/np.sqrt(np.pi))*qi**2
            for j in range(i):
                qj = charges[j]
                ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                rij = self.Rij[i,j]
                C += qi*qj*self.calc_Vscreen(ga_ij,rij)
        '''
        return C
                

    def min_Charges(self):
        from scipy.optimize import minimize
        q = self.q
        cons = ({'type': 'eq', 'fun': lambda x:  np.sum(x)-self.Qtot})
        res = minimize(self.fun_Etot, q, method='SLSQP', constraints=cons,
                       options={'maxiter': 500, 'ftol': 1e-08, 'iprint': 1, 'disp': True, 'eps': 1.4901161193847656e-08})
        #res = minimize(self.fun_Etot, q, method='COBYLA', constraints=cons,
        #               options={'maxiter': 500, 'ftol': 1e-08, 'iprint': 1, 'disp': True, 'eps': 1.4901161193847656e-08})

        #trust-constr
        #print(res.x)
        self.q = res.x

    def fun_Etot(self,q):
        #print("calculating electrostatic energy")
        #print(self.q)
        #self.E = 0.0
        E = 0.0
        for i in range(len(q)):
            #A = self.atoms[i]
            #para = self.par[A]
            qi = q[i]
            ga_ii = self.calc_gamma(self.atsize[i],self.atsize[i])

            E += self.eneg[i]*qi
            E += 0.5*(self.hard[i]+2.0*ga_ii/np.sqrt(np.pi))*qi**2

            for j in range(i):
                #B = self.atoms[j]
                #parb = self.par[B]
                qj = q[j]
                ga_ij = self.calc_gamma(self.atsize[i],self.atsize[j])
                rij = self.Rij[i,j]

                E += qi*qj*self.calc_Vscreen(ga_ij,rij)
        return E
