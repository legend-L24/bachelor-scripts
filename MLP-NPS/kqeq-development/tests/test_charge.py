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

mols = read("plusTesting.xyz@:",format='extxyz')
mols_train = mols[:18]
mols_test = mols[18:]
charge_train = [1 for _ in mols_train]
charge_test = [1 for _ in mols_test]

dipoles_prediction =  {
    0: [ 0.14456719,  0.99218286, -0.05241394,  0.15615882,  1.24292363, -1.11686282],
    1: [ 0.59702408, -0.25354497,  0.61470189,  0.12195736,  0.11741544, -1.00327248],
    2:[ 0.38976543,  0.00336875,  0.68947579,  0.08628019,  0.47116955, -0.65157503]}
charges_prediction = {
    0: [-0.35679261,  0.05844298, -0.15921459,  0.08334536, -0.18125668,  0.16437178,
 -0.12376036,  0.21919362,  0.23165292,  0.22174264,  0.1861191,   0.19176806,
  0.17317004,  0.29121773, -0.25389316, -0.0355681,   0.0396533,   0.14548146,
 -0.10786503,  0.29911908, -0.02993576,  0.24422648,  0.24420486,  0.23810355,
  0.21647332],
    1: [-1.21557456,  0.71843972, -0.48958019,  1.1258175,  -0.51682547,  0.48208362,
 -0.49120012,  0.32537955,  0.35884006,  0.34163464,  0.00223033, -0.0936657,
 -0.10497429,  0.5573949,   0.00795467, -0.15266642, -0.6726219,   1.85375798,
 -0.81045338,  0.87483693, -0.19104226,  0.13103099,  0.13099946,  0.10110456,
 -0.27290064],
    2:[-0.5748886,   0.15726611, -0.11905431,  0.25625286, -0.15125298,  0.14237518,
 -0.25118633,  0.21083764,  0.22759377,  0.21896364,  0.17159992,  0.10688473,
  0.09923572,  0.50537264, -0.05194895,  0.01696476, -0.12599785,  0.43166372,
 -0.15624799,  0.3129868,   0.04100041,  0.14049515,  0.14048553,  0.12428375,
  0.12631467]}
eneg_prediction = {
    0: [ 0.00423463,  0.00533132,  0.07836382, -0.00605065,  0.07867511,  0.01240064,
  0.01436924, -0.10228542, -0.10222835, -0.10219905, -0.10676287, -0.10554417,
 -0.10612304, -0.15641294,  0.00441689,  0.03887228,  0.04220398,  0.00568815,
  0.07875141,  0.01846107,  0.08094222, -0.10272494, -0.10272615, -0.10259187,
 -0.10582633],
    1: [ 0.02239793, -0.21876352, -0.1004299,  -0.30168808, -0.0999366,  -0.21964378,
 -0.1149958,  -0.16794671, -0.1674646,  -0.16741671, -0.17790025, -0.17526286,
 -0.17838306, -0.4627356,  -0.13059985, -0.13654761, -0.10642667, -0.42222448,
 -0.10084445, -0.2519643,  -0.10166964, -0.16975001, -0.16975429, -0.16906697,
 -0.17373542],
    2:[ 1.12609691e-01,  3.10031545e-02,  7.08405838e-02, -1.24935087e-02,
  6.76657763e-02,  3.07426978e-02,  2.73782985e-02, -4.06190112e-03,
 -4.87162672e-03, -4.28857571e-03, -7.16906962e-02, -4.38938681e-02,
 -4.61325580e-02, -3.01769138e-01,  3.29040855e-02,  3.28755191e-02,
  6.65023864e-02, -2.77895201e-02,  6.38498095e-02, -7.18165189e-05,
  5.77089792e-02, -8.66133508e-03, -8.67411859e-03, -8.65363779e-03,
 -6.94809916e-02]}

charges_calc = {
    0: [-0.35679261,  0.05844298, -0.15921459,  0.08334536, -0.18125668,  0.16437178,
 -0.12376036,  0.21919362,  0.23165292,  0.22174264,  0.1861191,   0.19176806,
  0.17317004,  0.29121773],
    1: [-1.21557456,  0.71843972, -0.48958019,  1.1258175,  -0.51682547,  0.48208362,
 -0.49120012,  0.32537955,  0.35884006,  0.34163464,  0.00223033, -0.0936657,
 -0.10497429,  0.5573949 ],
    2:[-0.5748886,   0.15726611, -0.11905431,  0.25625286, -0.15125298,  0.14237518,
 -0.25118633,  0.21083764,  0.22759377,  0.21896364,  0.17159992,  0.10688473,
  0.09923572,  0.50537264]
}
energy_calc = {
    0: 0.3812022332489132,
    1: -10.313293253804773,
    2:-0.020279388260854344}
dipole_calc = {
    0: [ 0.14456719,  0.99218286, -0.05241394],
    1: [ 0.59702408, -0.25354497,  0.61470189],
    2:[0.38976543, 0.00336875, 0.68947579]}
forces_calc ={
    0: [[-0.06492726, -0.36180127,  0.12939022],
 [ 0.47913872,  0.29187971, -0.16560937],
 [ 0.13220046,  0.01731984,  0.15614152],
 [-0.4298215,  0.20044814,  0.33640985],
 [-0.05833037,  0.12075381, -0.15240032],
 [ 0.1203822,  -0.14688209, -0.27085935],
 [ 0.22748758, -0.13732172,  0.47087232],
 [ 0.40138359,  0.01557817, -0.05371723],
 [-0.24103675,  0.06488177, -0.41212858],
 [-0.17412368,  0.19444621,  0.3652013 ],
 [-0.51643949, -0.11804449, -0.03094734],
 [-0.08239879, -0.66417075, -0.17192737],
 [ 0.42455032,  0.38201146, -0.11527866],
 [-0.21806503,  0.14090122, -0.08514697]]
,
    1: [[ -1.41087991, -11.55995036,   2.2680026],
 [ -9.84522481,   5.40125781,  -7.27233177],
 [  3.89595419,  -0.58723056,   4.05825286],
 [  6.05914366,  -3.22483885,  -4.6964836 ],
 [ -1.75614668,   2.02261521,  -4.68814918],
 [  4.07109286,  -3.29654106,   6.14178469],
 [  5.24362606,  -4.2343873,    3.06168127],
 [-12.57922401,   4.02451464,   0.40608908],
 [  7.34459552,   2.36087973,  10.28688859],
 [  5.36552905,  -0.25893565, -11.86570086],
 [  8.46565883,   1.43163422,  -0.97091829],
 [ -5.68658435,   0.50030633,   3.70872109],
 [ -4.3227137,    3.64517277,   3.86934196],
 [ -4.84482672,   3.77550339,  -4.30717822]],
    2:[[-1.20841367, -7.80609922,  1.2400541],
 [ 1.22267364,  2.23975427, -0.90028583],
 [ 0.25059872, -0.0993994,  0.16701492],
 [ 0.43912193, -0.1760221,  -0.35229386],
 [ 0.03973316,  0.02017917, -0.30507739],
 [ 0.74687057, -0.70156464, -1.61328715],
 [ 2.58208757, -1.82118036,  5.20706175],
 [ 2.78574331,  1.44721356, -0.44270314],
 [-1.1957754,   2.01543999, -2.4808611 ],
 [-0.8039944,   2.5319052,   2.00331394],
 [-1.19458402, -0.14287323,  0.02164371],
 [-1.22116291, -1.89824999,  0.10414638],
 [ 0.5876547,   2.05693994,  0.3298559 ],
 [-3.03055316,  2.333957,   -2.97858217]]}

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
                                training_set=mols_train,
                                training_system_charges=charge_train)
        for kQeq_target in ["charges", "dipole", ["dipole","charges"]]:
            if len(kQeq_target) == 2:
                wt = [1,0.5]
            else:
                wt = 1
            my_kqeq = kernel_qeq(Kernel=SOAP_Kernel,
                                    scale_atsize=atsize,
                                    radius_type=radtype)
            my_kqeq.train(targets=kQeq_target,lambda_reg=0.1,target_weights = wt)
            dipoles, charges, enegs = my_kqeq.predict(mols_test,predict_system_charges=charge_test)
            self.assertTrue(np.allclose(dipoles, np.array(dipoles_prediction[targ_count])))
            self.assertTrue(np.allclose(charges, np.array(charges_prediction[targ_count])))
            self.assertTrue(np.allclose(enegs, np.array(eneg_prediction[targ_count])))
            print(f"{kQeq_target} target is OK for simple SOAP predictions") 
            result = my_kqeq.calculate(mols_test[0], charge = 1)
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

