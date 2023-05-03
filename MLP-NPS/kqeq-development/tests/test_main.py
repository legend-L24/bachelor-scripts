import unittest
import sys
sys.path.append("..")
from ase.io import read
import numpy as np
from kqeq.kqeq import kernel_qeq
from kqeq.funct import get_dipoles
from kqeq.kernel import kernel
import matplotlib.pyplot as plt
#from HTMLTestRunner import HTMLTestRunner

mols = read("testcase.xyz@:",format='extxyz')
mols_train = mols[:18]
mols_test = mols[18:]

dipoles_prediction =  {
    0: [0.1894312,   0.46152806, -0.17109397,  0.09830022,  0.6487897,  -0.72382875],
    1: [0.50623716,  0.13551382, -0.23561542,  0.07082692,  0.45985141, -0.52700524],
    2:[0.37157295,  0.24587882, -0.12600048,  0.07499142,  0.49650463, -0.5560235 ]}
charges_prediction =  {
    0:[-0.38669309,  0.02061713, -0.26398943,  0.05111699, -0.33730083,  0.36503046,
               -0.45455385,  0.14617992,  0.14689721,  0.14392639,  0.10432088,  0.10385137,
               0.09685094,  0.26374592, -0.29777007, -0.18283351, -0.15659441,  0.1649488,
               -0.31032751,  0.50558116, -0.29939203,  0.15987315,  0.15987243,  0.15851995,
               0.09812202],
    1:[-1.24259312,  0.33477282, -0.33355448, -0.05462003, -0.34892114,  0.3647388,
              -0.82429861,  0.36380092,  0.36079069,  0.3665764,   0.12940403,  0.17274332,
              0.1689322,   0.5422282,  -1.10923868, -0.10776921, -0.34319092,  0.29795593,
              -0.36515661,  0.61406792, -0.31760053,  0.39904107,  0.39904413,  0.38925864,
              0.14358826],
    2:[-0.45588578,  0.09859634, -0.31387355,  0.10985514, -0.3516739,   0.33276772,
               -0.64807353,  0.15750468,  0.15645258,  0.15850779,  0.1226888,   0.10248428,
               0.09752308,  0.43312635, -0.34536938, -0.27512794, -0.26536075,  0.2114816,
               -0.29094763,  0.5713889,  -0.31594889,  0.19445126,  0.19445694,  0.18820461,
               0.13277128],}
eneg_prediction = {
    0: [ 0.03062744,  0.02626762,  0.10353519,  0.02644043,  0.10075072,  0.01177561,
                 0.06541253, -0.06334718, -0.06324099, -0.06323949, -0.0645706,  -0.06454343,
                 -0.06509997, -0.12236938,  0.0390454,   0.0859848,   0.08977085,  0.02795635,
                 0.09838393,  0.0113724,   0.09743359, -0.06366366, -0.06366434, -0.06353904,
                 -0.06376164],
    1: [ 0.07022971, -0.04506755,  0.05009422, -0.03454599,  0.04521525, -0.02314361,
                0.06193594, -0.1740102,  -0.17375145, -0.17374752, -0.17701642, -0.17700553,
                -0.17834917, -0.34465709,  0.09578554,  0.09747302,  0.09712889, -0.04010596,
                0.04047102, -0.05157481,  0.04025943, -0.174887,   -0.17488866, -0.17462231,
                -0.1751248 ],
    2:[ 0.02420579, 0.00256276,  0.08981505, -0.00337439,  0.08586379,  0.01834718,
                0.07964304, -0.07715502, -0.07732152, -0.07710668, -0.10364637, -0.09347728,
                -0.09495928, -0.24092962,  0.05002629,  0.12778929,  0.13498553,  0.02316388,
                0.08164629,  0.00741968,  0.09436029, -0.07902669, -0.07903218, -0.07880661,
                -0.1016026]}

charges_calc = {
    0:[-0.38669309,  0.02061713, -0.26398943,  0.05111699, -0.33730083,
               0.36503046, -0.45455385,  0.14617992,  0.14689721,  0.14392639,
               0.10432088,  0.10385137,  0.09685094,  0.26374592],
    1: [-1.24259312,  0.33477282, -0.33355448, -0.05462003, -0.34892114,
               0.3647388 , -0.82429861,  0.36380092,  0.36079069,  0.3665764 ,
               0.12940403,  0.17274332,  0.1689322 ,  0.5422282 ],
    2: [-0.45588578,  0.09859634, -0.31387355,  0.10985514, -0.3516739 ,
                0.33276772, -0.64807353,  0.15750468,  0.15645258,  0.15850779,
                0.1226888 ,  0.10248428,  0.09752308,  0.43312635]
}
energy_calc = {
    0: -2.3996511008863806,
    1:-8.879453867063244,
    2:-3.9106875565766606
}
dipole_calc = {
    0: [ 0.1894312 ,  0.46152806, -0.17109397],
    1: [ 0.50623716,  0.13551382, -0.23561542],
    2:[ 0.37157295,  0.24587882, -0.12600048]
}
forces_calc ={
    0:[[ 0.03890074,  0.52583954, -0.00231808],
               [ 0.13577025,  0.19346944,  0.34198442],
               [-0.23625035, -0.0948381,  -0.45461385],
               [ 0.12103633, -0.02004135, -0.1905929 ],
               [ 0.35486377, -0.27186083,  0.3742946 ],
               [-0.26986329,  0.24826227,  0.00672804],
               [ 0.30717929, -0.2311033,   0.2330489 ],
               [ 0.02430122, -0.1622905,   0.00684408],
               [-0.03950521, -0.16829063, -0.04254084],
               [-0.03639765, -0.12817533,  0.08972084],
               [-0.21077309, -0.03795993, -0.02516236],
               [-0.02948361, -0.21709739, -0.04719643],
               [ 0.14223996,  0.13844016, -0.03219476],
               [-0.30201836,  0.22564594, -0.25800165]],
    1:[[ 6.91308115e-01,  5.80641030e+00, -6.77732776e-01],
              [ 2.07932373e-01, -3.32698644e-01, -1.61801066e-01],
              [-1.08534326e+00,  1.80533958e-01, -5.66661425e-01],
              [ 5.66124263e-01, -1.84696356e-01, -4.20002156e-01],
              [ 4.91600979e-01, -4.45125908e-01,  9.77429576e-01],
              [-7.47887081e-01,  1.01130462e+00,  4.30486987e+00],
              [ 3.35918688e+00, -2.99482328e+00, -9.07413326e-01],
              [-6.02702148e+00, -7.76774064e-01,  5.74711614e-01],
              [ 2.99157310e+00, -1.66191974e+00,  5.00684589e+00],
              [ 2.09440976e+00, -2.85080398e+00, -4.82757539e+00],
              [-6.75621651e-03,  2.19663398e-01, -1.86818482e-03],
              [-9.26474967e-02, -6.46198256e-01, -1.63741666e-01],
              [ 4.22385865e-01,  4.64621020e-01, -1.03943937e-01],
              [-2.86486586e+00,  2.21050679e+00, -3.03311700e+00]],
    2:[[-0.10717554, -0.05699385,  0.04530509],
               [ 0.38238136,  0.67827484,  0.16593167],
               [-0.39184202, -0.1756886,  -0.79178974],
               [ 0.07793907, -0.02570089, -0.17223225],
               [ 0.79026368, -0.67205766,  0.43189598],
               [-0.4254653,   0.48914264,  1.72066712],
               [ 0.37071722, -0.32133029, -0.13615511],
               [ 1.03416047, -0.15109034, -0.08040919],
               [-0.55637431,  0.0209385,  -0.87825231],
               [-0.42293251,  0.25360932,  0.96601447],
               [-0.75348117, -0.11000034,  0.05434969],
               [-0.14206285, -0.95051661, -0.23591687],
               [ 0.61739023,  0.70158595, -0.14421294],
               [-0.47351833,  0.31982727, -0.94519562]]}

class TestKqeq(unittest.TestCase):
    
    def test_prediction(self):
        targ_count = 0
        desdict = {"nmax" : 2,
                   "lmax" : 2,
                   "rcut" : 2,
                   "sigma": 0.2,
                   "periodic": False}

        atsize = 1.0
        radtype = "qeq"

        SOAP_Kernel = kernel(Kernel='SOAP',
                             Descriptor='SOAP',
                             multi_SOAP=False,
                             descriptor_dict=desdict,
                             training_set=mols_train)
        for kQeq_target in ["charges", "dipole", ["dipole","charges"]]:
            if len(kQeq_target) == 2:
                wt = [1,0.5]
            else:
                wt = 1
            my_kqeq = kernel_qeq(Kernel=SOAP_Kernel,
                                    scale_atsize=atsize,
                                    radius_type=radtype)
            my_kqeq.train(targets=kQeq_target,lambda_reg=0.1,target_weights = wt)
            dipoles, charges, enegs = my_kqeq.predict(mols_test)
            self.assertTrue(np.allclose(dipoles, np.array(dipoles_prediction[targ_count])))
            self.assertTrue(np.allclose(charges, np.array(charges_prediction[targ_count])))
            self.assertTrue(np.allclose(enegs, np.array(eneg_prediction[targ_count])))
            print(f"{kQeq_target} target is OK for simple SOAP predictions") 
            result = my_kqeq.calculate(mols_test[0])
            self.assertTrue( np.allclose(result["charges"], np.array(charges_calc[targ_count])))
            self.assertTrue( np.allclose(result["energy"], np.array(energy_calc[targ_count])))
            self.assertTrue( np.allclose(result["dipole_vector"], np.array(dipole_calc[targ_count])))
            self.assertTrue( np.allclose(result["forces"], np.array(forces_calc[targ_count])))
            print(f"{kQeq_target} target is OK for simple SOAP calculations")
            targ_count += 1
        

        
if __name__ == '__main__':
    #runner=HTMLTestRunner(output="test_report")
    #unittest.main(testRunner=runner)
    unittest.main()

