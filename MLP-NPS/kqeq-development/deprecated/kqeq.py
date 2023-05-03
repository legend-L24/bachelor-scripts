import numpy as np
from ase.io import *
from ase.data import covalent_radii
from qeq import charge_eq
from sklearn.preprocessing import normalize
from dscribe.descriptors import SOAP


def block_diag(all_as,dim):
    A_bar = np.zeros((dim,dim))
    counti = 0
    countj = 0

    for i,a in enumerate(all_as):
        for ia in range(a.shape[0]):
            for ja in range(a.shape[1]):
                A_bar[ia+counti,ja+countj] = a[ia,ja]
        counti += a.shape[0]
        countj += a.shape[1]
    return A_bar

def get_descriptors(mols):
    ref_qs = []
    descriptors = []
    for mol in mols:
        ref_qs.extend(mol.get_initial_charges())
        descriptors.extend(mol.get_atomic_numbers())
    return ref_qs, descriptors
        
def build_A_bar(mols,scale_atsize=2.0):
    # Build A_bar
    all_as = []
    all_element = []
    ref_qs = []
    dim = 0
    for mol in mols:
        qe   = charge_eq(mol,scale_atsize=scale_atsize)
        dim += qe.nAt
        ref_qs.extend(mol.get_initial_charges())
        all_as.append(qe.get_Abar())
        all_element.extend(qe.atoms)

    A_bar = block_diag(all_as,dim)
    return A_bar, np.array(ref_qs), all_element

def build_A_bar_SOAP(mols,scale_atsize=3.0,species=["H", "C", "O", "N"],rcut=3.0,nmax=3,lmax=1):
    # Build A_bar
    all_as = []
    all_element = []
    ref_qs = []
    dim = 0
    soap = SOAP(
        species=species,
        periodic=False,
        rcut=rcut,
        nmax=nmax,
        lmax=lmax)
    for mol in mols:
        qe   = charge_eq(mol,scale_atsize=scale_atsize)
        dim += qe.nAt
        ref_qs.extend(mol.get_initial_charges())
        all_as.append(qe.get_Abar())
        #all_element.extend(qe.atoms)
        all_element.extend(normalize(soap.create(mol)))
    A_bar = block_diag(all_as,dim)
    return A_bar, np.array(ref_qs), all_element


def Kernel(A,B):
    if A==B:
        return 1.
    else:
        return 0.

def SOAP_Kernel(A,B):
    return np.dot(A,np.transpose(B))**2


def train_kQEq(mols,lams=np.logspace(-5,2,10),scale_atsize=2.0):
    A_bar, Q_ref, descriptors = build_A_bar(mols,scale_atsize=scale_atsize)
    
    K = np.zeros_like(A_bar)
    for i,ati in enumerate(descriptors):
        for j,atj in enumerate(descriptors):
            K[i,j] = Kernel(ati,atj)
    grid_qs = []
    lammin = 0
    minerror = 100.
    chinewmin = 0.0
    
    for il,lam in enumerate(lams):
        tmp = np.matmul(np.matmul(A_bar.T,A_bar),K)
        #print(tmp)
        tmp = np.linalg.inv(tmp+lam*np.eye(Q_ref.shape[0]))
        weights = np.matmul(np.matmul(tmp,A_bar.T),Q_ref)

        chinew = np.matmul(K,weights)
        qnew = np.matmul(A_bar,chinew)
        dQ = Q_ref-qnew
        mae = np.mean(np.abs(dQ))
        print(lam,mae)
        if mae < minerror:
            lammin = il
            minerror = mae
            chinewmin = chinew
        grid_qs.append(qnew)
    #print(weights)

    seen = []
    print('Fitted Electronegativities:')
    for i,At in enumerate(descriptors):
        if At not in seen:
            print(At,chinewmin[i])
            seen.append(At)
            
    return qnew,Q_ref,weights,descriptors
    #for i,q in enumerate(grid_qs[lammin]):
    #    print(q,Q_ref[i])
    
def Predict_Qs(mols,descriptors_train,weights,scale_atsize=2.0):
    A_bar, Q_ref, descriptors_predict = build_A_bar(mols,scale_atsize=scale_atsize)
    
    K = np.zeros((len(descriptors_predict),len(descriptors_train)))
    for i,ati in enumerate(descriptors_predict):
        for j,atj in enumerate(descriptors_train):
            K[i,j] = Kernel(ati,atj)
    
    chinew = np.matmul(K,weights)
    qnew = np.matmul(A_bar,chinew)
    dQ = Q_ref-qnew
    mae = np.mean(np.abs(dQ))
    print('MAE=',mae)
    return qnew,Q_ref



