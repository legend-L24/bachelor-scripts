import os
import sys
from multiprocessing import Process, Queue
import pandas as pd
import optuna
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time
from multiprocessing import Process, Pool
import numpy as np
import threading
import queue

sys.path.append('/home/jyb/.conda/envs/yfh/GSASII')

global exitFlag
exitFlag = 0

class RefineProject():
    def __init__(self,cif,exprm,wdir,csv):
        self.cif = cif 
        self.exprm = exprm 
        self.wdir = wdir 
        self.csv = csv
        #self.hash = str(hash(cif+exprm+wdir+csv))
    
    def putgpx(self,trial_number):
        import GSASIIscriptable as G2sc
        import shutil

        self.gpx = G2sc.G2Project(newgpx=os.path.join(self.wdir,'refine_trial_%d.gpx'%(trial_number)))
        self.hist = self.gpx.add_powder_histogram(self.csv,self.exprm)
        self.phase = self.gpx.add_phase(
            self.cif,
            phasename='trial',
            histograms=[self.hist]
        )
        self.hist.data['Instrument Parameters'][0]['I(L2)/I(L1)'] = [0.5, 0.5, 0]

        # Set to use iso
        for val in self.phase.data['Atoms']:
            val[9] = 'I'

    def refine_and_calc_Rwp(self, param_dict):
        self.gpx.do_refinements([param_dict])
        for hist in self.gpx.histograms():
            _, Rwp = hist.name, hist.get_wR()
        return Rwp

class OptunaProject():
    def __init__(self,cif,exprm,csv,max_time):
        self.cif = cif 
        self.exprm = exprm 
        self.csv = csv
        filename = csv.split('/')[-1]
        self.uid = filename
        self.max_time = max_time
        self.wdir = 'job_'+self.uid
        #print(self.cif+'   '+self.uid)
        os.system('echo %s >> Error_%s'%(self.cif,self.uid))
        

        if not os.path.exists(self.wdir):
            os.mkdir(self.wdir)
    
    # Get the lower and upper boundary of the spectrum
    def get_boundary(self):
        '''
        with open(self.csv,'r') as f:
            lines=f.readlines()
        
        lb = eval(lines[0].strip().split(',')[0])
        for i in range(len(lines)):
            if lines[-1-i].strip() != '':
                ub = eval(lines[-1].strip().split(',')[0])
                break

        print(lb,ub)
        '''
        lb=4.0
        ub=30.0
        return lb,ub

    def objective(self,trial):
        """
        
        Parameters
        ----------
        trial : optuna.trial object

        Returns
        -------
        Rwp : float
        
        """
        lb,ub = self.get_boundary()
        ### define search space ###
        # Limits (acute angle)
        limits_lb = trial.suggest_uniform('Limits lower bound', lb, lb*2/3+ub/3)
        limits_ub = trial.suggest_uniform('Limits upper bound', limits_lb + (ub-lb)/3, ub)
        limits_refine = trial.suggest_categorical('limits refine', [True, False])
        refdict0 = {'set': {'Limits': [limits_lb, limits_ub]}, 'refine': limits_refine}

        # Background
        background_type = trial.suggest_categorical(
            'Background type', ['chebyschev',
                                'cosine',
                                'Q^2 power series',
                                'Q^-2 power series',
                                'lin interpolate',
                                'inv interpolate',
                                'log interpolate'])
        no_coeffs = trial.suggest_int('Number of coefficietns', 1, 15 + 1)  # [1, 16)
        background_refine = trial.suggest_categorical('Background refine', [True, False])
        refdict0bg_h = {
            'set': {
                'Background': {
                    'type': background_type,
                    'no. coeffs': no_coeffs,
                    'refine': background_refine
                }
            }
        }
        # Instrument parameters
        instrument_parameters1_refine = []
        for p in ['Zero']:
            if trial.suggest_categorical('Instrument_parameters refine %s' % (p), [True, False]):
                instrument_parameters1_refine.append(p)
        refdict1_h = {'set': {'Cell': True, 'Instrument Parameters': instrument_parameters1_refine}}

        sample_parameters1_refine =[]
        for p in ['DisplaceX', 'DisplaceY']:
            if trial.suggest_categorical('Sample_parameters refine %s' % (p), [False]):
                sample_parameters1_refine.append(p)
        for p in ['Scale']:
            if trial.suggest_categorical('Sample_parameters refine %s' % (p), [True, False]):
                sample_parameters1_refine.append(p)
        refdict1_h2 = {"set": {'Sample Parameters':sample_parameters1_refine }}

        instrument_parameters2_refine = []
        for p in ['U', 'V', 'W', 'X', 'Y', 'SH/L']:
            if trial.suggest_categorical('Peakshape_parameters refine %s' % (p), [True, False]):
                instrument_parameters2_refine.append(p)
        refdict2_h = {'set': {'Instrument Parameters': instrument_parameters2_refine}}
        
        with open(self.cif,'r') as ciff:
            lines = ciff.readlines()
            for i in range(len(lines)):
                if '_atom_site_fract_z' in lines[i]:
                    starti = i
            atoms_dic = {}
            for i in range(len(lines)-starti-1):
                if lines[starti+i+1].rstrip() != '':
                    a = lines[starti+i+1].strip()
                    #with open('debug','w') as debugf:
                    #    debugf.write(a)
                    a = a.split()
                    if len(a) > 4:
                        a = a[0]
                    else:
                        break
                    for each in ['C','H','O','N']:
                        if each in a:
                            atoms_dic[a] = 'XU'
                        else:
                            atoms_dic[a] = 'XU'

            #print()
        refdict3_h = {'set': {'Atoms': {'all': 'XU'}}}
        #refdict3_h = {'set': {'Atoms': atoms_dic}}

        # Limits (wide angle)
        refdict_fin_h = {'set': {'Limits': [lb,ub]}, 'refine': True}

        # Evaluate
        refine_params_list = [refdict0,
                            refdict0bg_h,
                            refdict1_h,
                            refdict1_h2,
                            refdict2_h,
                            refdict3_h,
                            refdict_fin_h]

        def refine(trial_number,refine_params_list,q):
            ERROR_PENALTY = 1e9 
            gpxproject = RefineProject(self.cif,self.exprm,self.wdir,self.csv)
            gpxproject.putgpx(trial_number)
            Rwp = ERROR_PENALTY 
            for i in range(len(refine_params_list)):
                params = refine_params_list[i]
                Rwp_temp = gpxproject.refine_and_calc_Rwp(params) 
                if i==len(refine_params_list)-1:
                    Rwp = Rwp_temp 
                else: 
                    Rwp = ERROR_PENALTY
            # validate Uiso >= 0
            phase = gpxproject.gpx.phases()[0]
            u_iso_list = [atom.uiso for atom in phase.atoms()]
            if min(u_iso_list) < 0:
                # Uiso < 0
                Rwp = ERROR_PENALTY
            with open('temp_'+self.uid,'w') as f:
                f.write('%.10f'%Rwp)

        def evaluate(trial_number, refine_params_list, q):
            p1 = Process(target=refine,args=(trial_number,refine_params_list,q,))
            p1.start()
            p1.join(self.max_time)
            if p1.is_alive():
                p1.terminate()
                with open('temp_'+self.uid,'w') as f:
                    f.write('1000000000')
                    
            with open('temp_'+self.uid,'r') as f:
                temperror = eval(f.readline())
                q.put(temperror)
                os.system('echo %f >> Error_%s'%(temperror,self.uid))

        self.q = Queue()
        self.p = Process(target=evaluate, args=(trial.number, refine_params_list, self.q))
        self.p.start()
        Rwp = self.q.get()
        self.p.join()

        return Rwp

    def create_and_optimize(self,random_seed,ntrials):
        self.study = optuna.create_study(study_name=self.cif + '_seed%s' % (random_seed),
                                         sampler=optuna.samplers.TPESampler(n_startup_trials=20, seed=random_seed))
        self.study.optimize(self.objective, n_trials=ntrials, n_jobs=1)

class myThread(threading.Thread):
    def __init__(self,name,q):
        threading.Thread.__init__(self)
        self.name = str(name) 
        self.q = q


    def run(self):
        process_data(self.q)

def process_data(q):
    while not exitFlag:
        if not q.empty():
            data = q.get()
            optproject = OptunaProject(data[0],data[1],data[2],data[3])
            optproject.create_and_optimize(data[4],data[5])
        #optproject = OptunaProject(cif,exprm,csv,max_time)
        #optproject.create_and_optimize(rand_seed,ntrials)


class MultiRefine():
    
    def __init__(self,cifpath,exprm,csv,rand_seed=0,ntrials=50,max_time=150,nthreads=3):
        self.cifpath = cifpath
        self.exprm = exprm
        self.csv = csv
        self.cif_list = os.listdir(self.cifpath)
        self.rand_seed = rand_seed 
        self.ntrials = ntrials
        self.max_time = max_time
        self.nthreads = nthreads
 
    def start(self):
        workQueue = queue.Queue(len(self.cif_list))
        # Create threads 
        threads_list = []
        for i in range(self.nthreads):
            thread = myThread(i,workQueue)
            thread.start()
            threads_list.append(thread)

        # Fill in queue 
        for ciffile in self.cif_list:
            workQueue.put([os.path.join(self.cifpath,ciffile),self.exprm,self.csv,self.max_time,self.rand_seed
            ,self.ntrials])

        # Wait for queue to be empty 
        while not workQueue.empty():
            time.sleep(10) # Check every 10 seconds
            pass 

        # Stop all threads after the queue is empty
        exitFlag = 1 
        for t in threads_list: 
            t.join() 


