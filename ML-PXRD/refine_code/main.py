from functions_refine import MultiRefine 

kwargs = {
    'cifpath'   : '/home/jyb/yifanhou/examples/1/cif',        # path to all cifs
    'exprm'     : '/home/jyb/yifanhou/examples/1/inst_xry.prm',    # instrument parameters
    'csv'       : '/home/jyb/yifanhou/examples/1/spec.csv',     # PXRD spectrum
    'rand_seed' : 1,                                                    # random seed
    'ntrials'   : 200,                                                   # Number of trials in each optuna work
    'max_time'  : 120,                                                  # Maximal time for BBO refinement
    'nthreads'  : 1                                                     # Number of threads
}

project = MultiRefine(**kwargs)
project.start()
