
#############################################################################
# uti_math_opt module: optimization utilities / functions
# v 0.03
# Authors: An He, O.C.
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from array import *
from math import *
from copy import *
import os
import errno
import datetime
#import tempfile
#import shutil
import numpy as np
#import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import shelve, json, logging
import mpi
import json

#****************************************************************************

def optimize(_f, _x, _x_lim, _meth, _fn=None, _opt=None, _jac=None, _hess=None, _eps=1.e-06,  _aux=None):
    """
    Wrapper of minimization methods that can be used with (limited-accuracy) forward-simulation codes such as SRW
    :param _f: the objective function to be minimized.
    :param _x: an array / list of instant values of optimization variables
    :param _x_lim: limiting (i.e. min. amd max. possible) values of parameters to be optimized
    :param _meth: optimization method number

    :return OptimizeResult. The optimization result represented as a OptimizeResult object. 
            Important attributes are: x the solution array, success a Boolean flag indicating if the optimizer exited successfully and
            message which describes the cause of the termination.
    """
    if _meth == 0:   # scanning
        opt = {'grid': [2,2,2,2], 'sqc': 0, 'samp': 4}
        if(_opt is not None): opt.update(_opt)
        return optimize_grid_scan(_f, _opt, _fn, (_x_lim, _aux))

    elif _meth == 1: # polynomial fitting #AH02202019
        print('Start to polynomial fit')
        opt = {'deg': 3, 'samp':4}
        if(_opt is not None): opt.update(_opt)
        if len(_x_lim) > 1: raise ValueError('Only one parameter is allowed to do polynomial fitting.') 
        if _opt['deg'] > _opt['samp']: raise ValueError('Degree of the fitting polynomial can not be larger than the number of sample points.')
        return optimize_polynomial_fit(_f, _opt, args=(_x_lim, _aux))

    elif _meth == 2: #'downhill simplex': #AH01242019
        opt = {'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': None, 'maxfun': None, 'full_output': 0, 'disp': 1, 'retall': 0, 'initial_simplex': None}
        if(_opt is not None): opt.update(_opt)
        return scipy.optimize.fmin(_f, _x, args=(_x_lim, _aux), xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0)

    elif _meth == 3: #'Powell':
        opt = {'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': None, 'maxfev': None, 'disp': True, 'direc': None, 'return_all': True}
        if(_opt is not None): opt.update(_opt)
        return scipy.optimize.minimize(_f, _x, args=(_x_lim, _aux), method='Powell', tol=None, callback=None, options=opt)

    elif _meth == 4: #'L-BFGS-B':
        #res = L_BFGS_B(_f, _x, jac = _jac, bound = _bound, epsilon = _epsilon)
        #jac = kwargs['jac']
        #bound = kwargs['bound']
        #epsilon = kwargs['epsilon']
        opt = {'disp': True, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-05,
               'maxfun': 15000, 'maxiter': 15000, 'iprint': 99, 'maxls': 20}
        if(_opt is not None): opt.update(_opt)
        return scipy.optimize.minimize(_f, _x, args=(_x_lim, _eps, _f), method='L-BFGS-B', jac=_jac, bounds=_bnds, tol=None, callback=None, options=opt)

    elif _meth == 5: #'trust-exact':
        #res = trust_exact(_f, _x, jac = _jac, hess = _hess, bound = _bound, epsilon = _epsilon)
        #if scipy.__version__ != '1.1.0' : raise Exception ('Warning: SciPy v1.1.0 is required for the optimization method "trust_exact"!')
        #jac = kwargs['jac']
        #hess = kwargs['hess']
        #bound = kwargs['bound']
        #epsilon = kwargs['epsilon']
        opt = {'gtol': 1e-8, 'disp': True}
        if(_opt is not None): opt.update(_opt)
        return scipy.optimize.minimize(_f, _x, args=(_x_lim, _jac, _eps, _f), method='trust-exact', jac=_jac, hess=_hess, tol=None, callback=None, options=opt)

    else: return None

#****************************************************************************
def optimize_grid_scan(f, opt, fn, args):
    """The function does least squares polynomial fit, and then get the extrema of the polynome
    :param opt['grid']: an array, the grid size of each dimension
    :param opt['samp']: number of scanning points
    :param opt['sqc']: perform sequential calculation if sqc=0, parallel calculation if sqc=1
    """
    fcalls, func = wrap_function(f, args)
    xlim = args[0]
    grid = opt['grid']
    sqc = opt['sqc']
    num = opt['samp']

    x = [np.linspace(xlim[i][0], xlim[i][1], num*grid[i]) for i in range(len(grid))]
    
    import itertools
    pop = list(itertools.product(*x)) # to produce the list of the grid according the dimention and the size. 
    if sqc == 1:  # do the sequential calculation to obtain the cost function of each points
        pop = [{'index': i, 'x': pi, 'costf': func(np.array(pi))} for i,pi in enumerate(pop)]
        #print('pop',pop)
        #with open('data_optimize_grid_scan.txt', 'w') as f:
        with open(fn, 'w') as f:
            json.dump(pop, f, indent=2) # save the sample points to data file.
    elif sqc == 0: # do the parallel calculation to obtain the cost function of each points. and the data is save through mpi.run()
        pop = [np.array(pi) for pi in pop]
        dat = mpi.run(np.array(pop), func)
    return
#****************************************************************************
def wrap_function(function, args): #AH02202019
    ncalls = [0]
    if function is None:
        return ncalls, None

    def function_wrapper(*wrapper_args):
       ncalls[0] += 1
       return function(*(wrapper_args + args))

    return ncalls, function_wrapper

#****************************************************************************
def optimize_polynomial_fit(f, opt, args): #AH02202019
    """The function does least squares polynomial fit, and then get the extrema of the polynome
    :param opt['deg']: degree of fitting polynomial
    :param opt['samp']: number of sample points
    """
    fcalls, func = wrap_function(f, args)
    # initialize x and f
    xlim = args[0]
    left, right = xlim[0][0], xlim[0][1]
    x = np.linspace(left, right, opt['samp'])
    f = [func(np.array([xi])) for xi in x]   
    #np.savetxt('data_optimize_polynomial_fit.txt',[x,f])   # save the sample points to data file.
    z = np.polyfit(x, f, opt['deg']-1) # Least squares polynomial fit
    p = np.poly1d(z)  # get the polynome of z
    d = p.deriv()    # return a derivative of this polynomial
    mini = np.poly1d(d).r # return the polynome and its extrema
    x1 = np.linspace(left, right, 100)
    plt.plot(x,f,'.', x1,p(x1),'--')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('f')
    plt.savefig('poly_fit_plot.png')
    print('save and check the result:')
    import cmath
    res = [[m.real, func(np.array([m.real])), p(m)] for m in mini if m.imag == 0] # res[0]:x_min, res[1]:cost_fun(x_min), res[2]:polynominal(x_min)
    if len(res)==0:
        raise Exception('no real solution for this polynomial fitting, please change the number of fitting points or the degree of polynomial fitting. Try again!')
    else:
        print('polynomial fitting result:')
        for i in res:
            print('parameters = {}, cost_func(x_min) = {}, p(x_min)'.format(i[0], i[1]))
    return res

#****************************************************************************
def check_limits(_x, _lim):
    """The function checks if variables _x are between the limits specified by _lim
    :param _x: an array / list of instant values of optimization variables
    :param _lim: an array / list of pairs defining acceptable lower and upper limits of optimization variables
    :return: True for _x between limits _lim, False otherwise
    """

    len_x = len(_x)
    len_lim = len(_lim)
    if(len_x > len_lim): raise Exception("Inconsistent numbers of optimization parameters and limits limits")
    
    for i in range(len_x):
        #if _x[i] < _lim[i,0] or _x[i] > _lim[i,1]: return False
        curLim = _lim[i]
        if _x[i] < curLim[0] or _x[i] > curLim[1]: return False
    return True

#****************************************************************************
def norm_weights(_w):
    """Normalize weights for the optimization.
    :param _w: input non-normalized weights
    :return: normalized weights
    """
    nWeights = len(_w)
    sumWeights = 0
    for i in range(nWeights): sumWeights += _w[i]

    if(sumWeights == 0): return 0

    invSum = 1./sumWeights
    for i in range(nWeights): _w[i] *= invSum
    
    return _w

#****************************************************************************
def status_str(_x_nm, _x_inst, _cost_inst, _num_frm='{:04.6g}', _stat=None):
#def status_str(_x_nm, _x_inst, _cost_inst, _num_frm='{:04.6g}', _x_best=None, _cost_best=None, _call_num=None):
    """Generate instant optimization status string.
    :param _x_nm: list of names of optimization parameters
    :param _x_inst: list of instant values of optimization parameters + the corresponding cost value
    :param _cost_inst: instant cost value
    :param _num_frm: format to use for representing all numbers
    :param _x_best: list of best (so far) values of optimization parameters + the corresponding cost value
    :param _cost_best: best (so far) cost value
    :param _call_num: cost function call number
    :return: formatted string describing optimization status
    """

    nPar = len(_x_nm)

    resStr = ''
    if(_stat is not None):
        if(len(_stat) < (nPar + 2)):
            raise Exception("Inconsistent number of optimization parameters and status list length")
        
        resStr += '   Cost function call #{:d}\n'.format(_stat[nPar + 1] + 1)

    resStr += '   Instant values of optimization parameters:\n'
    resStr += '  '
    for i in range(nPar): resStr += ' ' + _x_nm[i] + '=' + _num_frm.format(_x_inst[i])
    resStr += '  cost: ' + _num_frm.format(_cost_inst)# + '\n'

    if((_stat is not None) and (_stat[nPar + 1] > 0)):
        #if((_x_best is not None) or (_cost_best is not None)):
        resStr += '\n   Best (so far) values of optimization parameters:\n'
        #if(_x_best is not None):
        resStr += '  '
        for i in range(nPar): resStr += ' ' + _x_nm[i] + '=' + _num_frm.format(_stat[i])
        #if(_cost_best is not None):
        resStr += '  cost: ' + _num_frm.format(_stat[nPar])

    return resStr

#****************************************************************************
def status_update(_stat, _x, _cost):
    """Update optimization status array (i.e. current best values of optimization parameters, best cost, function call number)
    :param _stat: status list/array, starting from best values of parameters, followed by cost, followed by call number
    :return: updated status parameters
    """

    nPar = len(_x)
    if(len(_stat) < (nPar + 2)):
        raise Exception("Inconsistent number of optimization parameters and status list length")

    prevCallNum = _stat[nPar + 1]
    prevCost = _stat[nPar]
    if((_cost < prevCost) or (prevCallNum <= 0)):
        for i in range(nPar): _stat[i] = _x[i]
        _stat[nPar] = _cost
    _stat[nPar + 1] += 1
    
    return _stat

#****************************************************************************
def log_update(_fpath, _str):
    """Updates log file(s) of an optimization process
    """

    timestamp = '[{:%Y-%m-%d %H:%M:%S}]:'.format(datetime.datetime.now())

    #This leaves only one (last) record in the file:
    #tmp_file = tempfile.NamedTemporaryFile(delete=False, mode='a', dir=os.path.dirname(_fpath))
    #tmp_file.write(timestamp + _str + '\n')
    #tmp_file.close()
    #shutil.move(tmp_file.name, _fpath)
    
    with open(_fpath, 'a') as f: f.write(timestamp + _str + '\n')

##    # Save JSON file for Sirepo:
##    status = {
##        'timestamp': timestamp,
##        'particle_number': particle_number,
##        'total_num_of_particles': total_num_of_particles,
##        'progress': progress,
##        'status': status,
##    }
##    status_json_file = '{}.json'.format(filename)
##    #with open(status_json_file, 'w') as f:
##    #    json.dump(status, f, indent=4, separators=(',', ': '), sort_keys=True)
##    #MR05122016: A temp file is created on the same filesystem as the status file to ensure the atomic rename operation:  
##    # See the following discussions:  
##    # - https://github.com/radiasoft/sirepo/issues/555  
##    # - https://github.com/radiasoft/SRW-light/commit/987547f4da3079626bef92aad2f9ca3bade84ada  
##    # - http://stackoverflow.com/a/3716361/4143531  
##    tmp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=os.path.dirname(status_json_file))
##    json.dump(status, tmp_file, indent=4, separators=(',', ': '), sort_keys=True)
##    tmp_file.close()
##    shutil.move(tmp_file.name, status_json_file)    
    
#****************************************************************************
def log_init(_dir, _str=None):
    """Creates log file(s) for auxiliary information generated during optimization
    :param _dir: log directory
    :param _str: string to be saved in the log file
    :return: path to the log file created
    """
    #This is mainly copied from srwlib.srwl_uti_save_stat_wfr_emit_prop_multi_e_init by MR
    
    log_dir = os.path.abspath(_dir) #os.path.abspath('__srwl_logs__')  
    #if rank == 0:  
    try:
        os.mkdir(log_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(log_dir): pass
        else: raise
        
    timestamp = '{:%Y-%m-%d_%H-%M-%S}'.format(datetime.datetime.now())
    log_file = 'uti_math_opt_minimize_{}.log'.format(timestamp)
    log_path = os.path.join(log_dir, log_file)

    strInit = timestamp = '[{:%Y-%m-%d %H:%M:%S}]:   Optimization started\n'.format(datetime.datetime.now())
    with open(log_path, 'w') as f: f.write(strInit)

    if(_str is not None): log_update(log_path, _str)
    #srwl_uti_save_stat_wfr_emit_prop_multi_e(0, total_num_of_particles, filename=log_path, cores=num_of_proc, particles_per_iteration=num_part_avg_proc)  

    return log_path



#****************************************************************************
#****************************************************************************
#****************************************************************************
def cal_real_num(_n_part_tot, _n_part_avg_proc, _nProc):
    """
    For partially coherent, to caculate the real number of macro-electrons to be calculated (for parallel calculation)
    :param _n_part_tot: total number of macro-electrons to be set
    :param _n_par_avg_proc: number of macro-electrons / wavefront to average on worker processes before sending data to master
    :param _nProc: number of cores
    :return: the real number of macro-electrons to be calculated
    """
    nPartPerProc_1 =int(round( _n_part_tot//(_nProc - 1)))  #Number of sending acts made by each worker process
    nSentPerProc = int(round(nPartPerProc_1//_n_part_avg_proc)) #Number of sending acts made by each worker process
    nPartPerProc_2 = _n_part_avg_proc*nSentPerProc  #Number of electrons treated by each worker process

    return nPartPerProc_2*(_nProc-1)

#****************************************************************************
def check_bound(_x, _y):
    """The function checks if variables _x are between the bounds _y.
    :param _x: an array of _x values.
    :param _y: an array with np.shape(_y) = (len(_x),2).
    :return: True for _x between _y, False for _x out of bounds _y.
    """
    res = True
    for i in range(len(_x)):
        if _x[i] < _y[i,0] or _x[i] > _y[i,1]: res = False

    return res

#****************************************************************************
def J(_x, *args):
    """ To compute the gradient vector by call "Finite-difference approximation of the gradient of a scalar function" scipy.optimize.approx_fprime() .
    :param _x: array_like. The coordinate vector at which to determine the gradient of f.
    :param args[-1]: callable function f. The function of which to determine the gradient (partial derivatives). Should take _x as first argument, \
                     other arguments to f can be supplied in *args. Should return a scalar, the value of the function at _x.
    :param args[-2]: array_like. Increment to _x to use for determining the function gradient. If a scalar, uses the same finite difference delta \
                     for all partial derivatives. If an array, should contain one value per element of _x.
    :param *args[:-2]: optional. Any other arguments that are to be passed to f.
    :return: ndarray. The partial derivatives of f to _x.
    """

    f = args[-1]
    epsilon = args[-2]

    return scipy.optimize.approx_fprime(_x, f, epsilon, *args[:-2])

#****************************************************************************

def Hess(_x, *args):
    """ To compute the Hessian matrix.
    :param _x: array_like. The coordinate vector at which to determine the gradient of f.
    :param _J: callable function fprime.
    :param args[-1]: callable function f. The function of which to determine the gradient (partial derivatives). Should take _x as first argument, \
                     other arguments to f can be supplied in *args. Should return a scalar, the value of the function at _x.
    :param args[-2]: array_like. Increment to _x to use for determining the function gradient. If a scalar, uses the same finite difference delta \
                     for all partial derivatives. If an array, should contain one value per element of _x.
    :param *args[:-2]: optional. Any other arguments that are to be passed to f.
    :return: Hessian Matrix (n, n).
    """
    _J = args[-3]
    g0 = np.zeros_like(_x, dtype='f')
    g0 = _J(_x, *args)
    gp = np.zeros([len(_x),len(_x)], dtype='f')
    x1 = 1.0 * np.copy(_x)
    epsilon = args[-2]
    g1 = np.zeros_like(_x, dtype='f')
    for i in range(len(_x)):
        x1[i] += epsilon[i]
        g1 = _J(x1, *args)
        print('g1:',g1)
        gp[i] = (g1-g0)/epsilon[i]
        x1[i] -= epsilon[i]

    return gp

#****************************************************************************
def Powell(_f, _x, **kwargs):
    """
    Method Powell is a modification of Powell’s method which is a conjugate direction method. 
    It performs sequential one-dimensional minimizations along each vector of the directions set (direc field in options and info), 
    which is updated at each iteration of the main minimization loop. The function need not be differentiable, and no derivatives are taken.
    :param _f: the objective function to be minimized.
               _f(x, *args) -> float
               where x is an 1-D array with shape (n,) and args is a tuple of the fixed parameters needed to completely specify the function.
    :param _x: Initial guess. ndarray, shape (n,)
               Array of real elements of size (n,), where ‘n’ is the number of independent variables.
    :param *kwargs: tuple, optional
                   Extra arguments passed to the objective function and its derivatives (fun, jac and hess functions).
    :return OptimizeResult. The optimization result represented as a OptimizeResult object. 
            Important attributes are: x the solution array, success a Boolean flag indicating if the optimizer exited successfully and
            message which describes the cause of the termination.
    """

    bound = kwargs['bound']

    return scipy.optimize.minimize(_f, _x, args=(bound,), method='Powell', tol=None, callback=None, options={'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': None, 'maxfev': None, 'disp': True, 'direc': None, 'return_all': True})

#****************************************************************************
def L_BFGS_B(_f, _x, **kwargs):
    """
    Method L-BFGS-B uses the L-BFGS-B algorithm for bound constrained minimization
    :param _f: the objective function to be minimized.
               _f(x, *args) -> float
               where x is an 1-D array with shape (n,) and args is a tuple of the fixed parameters needed to completely specify the function.
    :param _x: Initial guess. ndarray, shape (n,)
               Array of real elements of size (n,), where ‘n’ is the number of independent variables.
    :param *kwargs: tuple, optional
                   Extra arguments passed to the objective function and its derivatives (fun, jac and hess functions).
    :return OptimizeResult. The optimization result represented as a OptimizeResult object. 
            Important attributes are: x the solution array, success a Boolean flag indicating if the optimizer exited successfully and
            message which describes the cause of the termination.
    """

    jac = kwargs['jac']
    bound = kwargs['bound']
    epsilon = kwargs['epsilon']

    return scipy.optimize.minimize(_f, _x, args=(bound, epsilon, _f), method='L-BFGS-B', jac=jac, bounds=bound, tol=None, callback=None, options={'disp': True, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-05, 'maxfun': 15000, 'maxiter': 15000, 'iprint': 99, 'maxls': 20})
#    return scipy.optimize.minimize(_f, _x, args=(bound, epsilon), method='L-BFGS-B', jac=jac, bounds=bound, tol=None, callback=None, options={'disp': True, 'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-05, 'maxfun': 15000, 'maxiter': 15000, 'iprint': 99, 'maxls': 20})
#****************************************************************************
def trust_exact(_f, _x, **kwargs):

    """
    Method trust-exact is a trust-region method for unconstrained minimization in which quadratic subproblems are solved almost exactly. 
    This algorithm requires the gradient and the Hessian (which is not required to be positive definite). It is, in many situations, 
    the Newton method to converge in fewer iteraction and the most recommended for small and medium-size problems.
    :param _f: the objective function to be minimized.
               _f(x, *args) -> float
               where x is an 1-D array with shape (n,) and args is a tuple of the fixed parameters needed to completely specify the function.
    :param _x: Initial guess. ndarray, shape (n,)
               Array of real elements of size (n,), where ‘n’ is the number of independent variables.
    :param *kwargs: tuple, optional
                   Extra arguments passed to the objective function and its derivatives (fun, jac and hess functions).
    :return OptimizeResult. The optimization result represented as a OptimizeResult object. 
            Important attributes are: x the solution array, success a Boolean flag indicating if the optimizer exited successfully and
            message which describes the cause of the termination.
    """
    if scipy.__version__ != '1.1.0' : raise Exception ('Warning: SciPy v1.1.0 is required for the optimization method "trust_exact"!')
    jac = kwargs['jac']
    hess = kwargs['hess']
    bound = kwargs['bound']
    epsilon = kwargs['epsilon']

    return scipy.optimize.minimize(_f, _x, args=(bound, jac, epsilon, _f), method='trust-exact', jac=jac, hess=hess, tol=None, callback=None, options={'gtol': 1e-8, 'disp': True})

#****************************************************************************
def Jac(_fun, _x, *args): # didn't use it any more
    print('call Jac')
    print('*args=',args)
    f0 = _fun(_x, *args)
    fp = np.zeros_like(_x, dtype='d')
    x1 = 1.0 * np.copy(_x)
    dx = args[1]
    for i in range(len(_x)):
        x1[i] += dx[i]
        print('jjjjjac_x',x)
        f1 = _fun(x1, *args)
        fp[i] = (f1-f0)/dx[i]
        x1[i] -= dx[i]
    return fp

#****************************************************************************
def srwl_optimize(_f, _x, _method, _bound, _jac, _hess, _epsilon):
    """
    Method L-BFGS-B uses the L-BFGS-B algorithm for bound constrained minimization
    :param _method: the optimize method 

    :return OptimizeResult. The optimization result represented as a OptimizeResult object. 
            Important attributes are: x the solution array, success a Boolean flag indicating if the optimizer exited successfully and
            message which describes the cause of the termination.
    """
    if _method == 'Powell':
        res = Powell(_f, _x, bound = _bound)
    elif _method == 'L-BFGS-B':
        res = L_BFGS_B(_f, _x, jac = _jac, bound = _bound, epsilon = _epsilon)
    elif _method == 'trust-exact':
        res = trust_exact(_f, _x, jac = _jac, hess = _hess, bound = _bound, epsilon = _epsilon)
    return res
#****************************************************************************


#****************************************************************************
def sv_shelve(fn, *args):
    """To save the optimized results.
    :param fn: the file name to save the optimized result, ex: Optimizeresult.shelve.
    :param *args: any information which needed to pickle to the file fn.

    Help on module shelve:

    NAME
        shelve - Manage shelves of pickled objects.

    DESCRIPTION
    A "shelf" is a persistent, dictionary-like object.  The difference
    with dbm databases is that the values (not the keys!) in a shelf can
    be essentially arbitrary Python objects -- anything that the "pickle"
    module can handle.  This includes most class instances, recursive data
    types, and objects containing lots of shared sub-objects.  The keys
    are ordinary strings.
    
    To summarize the interface (key is a string, data is an arbitrary
    object):
    
            import shelve
            d = shelve.open(filename) # open, with (g)dbm filename -- no suffix
    
            d[key] = data   # store data at key (overwrites old data if
                            # using an existing key)
            data = d[key]   # retrieve a COPY of the data at key (raise
                            # KeyError if no such key) -- NOTE that this
                            # access returns a *copy* of the entry!
            del d[key]      # delete data stored at key (raises KeyError
                            # if no such key)
            flag = key in d # true if the key exists
            list = d.keys() # a list of all existing keys (slow!)
    
            d.close()       # close it
    """
    shvf = shelve.open(fn, 'c')
    shvf["res"] = args[0]
    shvf.close()
    return
#****************************************************************************
def open_shelve(fn, *args):
    """To read the shelve files.
    :param fn: the file name to be opened, ex: Optimizeresult.shelve.
    :param *args: any information which needed to pickle to the file fn.
    """
    res = shelve.open(fn, 'r')
    res.close() 
    return res
'n_fun	CostFunction	Targetparameters	Optimizeparameters'
#****************************************************************************
def save_opt_header_to_log(_target, _variable, *args):
    """save the header of the optimize information to the .log file.
    :param _target: list of the name of the optimization term.
    :param _variable: list of the name of the optimization variable.
    :param *args: any information which will be wrote to .log file.
    """
    logging.info('OPTIMIZATIOB RECORD')
    logging.info('optimization function:fun = xFWHMwn*((xFWHM - xFWHMtarg)/xFWHMnom)**2 + yFWHMwn*((yFWHM - yFWHMtarg)/yFWHMnom)**2\n                                                  + xCOHLwn*((xCOHL - xCOHLtarg)/xCOHLnom)**2 + yCOHLwn*((yCOHL - yCOHLtarg)/yCOHLnom)**2')
    logging.info('optimization term:	{}'.format(_target))
    logging.info('optimization variable:	{}'.format(_variable))
    logging.info('n_fun	= number of function called')
    logging.info('fun	= value of optimization function')
    logging.info('term	= optimization term')
    logging.info('x	= optimization variable')
    logging.info('MIN	= current minimum')
    logging.info('')
    header2_1 = ''.join(['{:^12}'.format(i) for i in _target])
    header2_2 = ''.join(['{:^12}'.format(i) for i in _variable])

    logging.info('{0:^4}{1:^12}{2:^52}{3:^74}{4:^3}'.format('n_fun','fun', 'term', 'x', 'MIN'))
    logging.info('{0:^4}{1:^12}{2:^52}{3:^74}'.format('','', header2_1, header2_2))
    logging.info(143*'-')
#****************************************************************************
def save_opt_to_json(_filename, _res, **kwargs):
    """save the initial and final optimized information to the .json file.
    :param _filename: the name without extension used to save log/status files.
    :param _res: optimized result from scipy.minimize.
    :param *kwargs: initial optimize information.
    """
    timestamp = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    x0 = kwargs['x0']	# initial variables
    f0 = kwargs['f0']	# initial function
    method = kwargs['method']	# optimize methods

    # Save JSON file for Sirepo:

    status = {
        'timestamp': timestamp,
        'x_ini': x0.tolist(),
        'f_ini': f0,
        'optimizer': method,
#        'success': _res.success,
#        'status': _res.status,
#        'message': _res.message,
#        'nfun': _res.nfev,
        'f_opt': _res.fun,
        'x_opt': _res.x.tolist(),
        'nit': _res.nit,
    }
    status_json_file = '{}.json'.format(_filename)
    f = open(status_json_file,'w')
    json.dump(status, f, indent=4, separators=(',', ': '), sort_keys=True)
    f.close()

#****************************************************************************
def save_opt_path():   
    """Initialize parameters for the SRW status files and generate the files.  
  
    """  
    log_dir = os.path.abspath('__srwl_opt_logs__')   
    try:  
        os.mkdir(log_dir)  
    except OSError as exc:  
        if exc.errno == errno.EEXIST and os.path.isdir(log_dir):
            pass  
        else:
            raise 
    timestamp = '{:%Y-%m-%d_%H-%M-%S}'.format(datetime.datetime.now())  
    log_file = 'srwl_optimize_{}'.format(timestamp)  
    log_path = os.path.join(log_dir, log_file)  

    return log_path
#****************************************************************************

#****************************************************************************
def fwhm(x, y, shift=0.5, return_as_dict=False):  # MR21032017
    """The function searches x-values (roots) where y=0 (after normalization to values between 0 and 1 and shifting the
    values down by 0.5 (default value)) based on linear interpolation, and calculates full width at half maximum (FWHM).

    :param x: an array of x values.
    :param y: an array of y values.
    :param shift: an optional shift to be used in the process of normalization (between 0 and 1).
    :param return_as_dict: if to return a dict with 'fwhm' and 'x_range'
    :return: a value of the FWHM or dictionary consisting of 'fwhm' and 'x_range'
    """

    def is_positive(num):
        return True if num > 0 else False

    # Normalize values first:
    #y = (y - min(y)) / (max(y) - min(y)) - shift  # roots are at Y=0

    #OC18112017 (making it work with standard Python lists / arrays)
    minY = min(y)
    maxY = max(y)
    if(maxY == minY): raise Exception('FWHM can not be calculated')
    mult = 1./(maxY - minY)
    lenY = len(y)
    for i in range(lenY): y[i] = (y[i] - minY)*mult - shift

    positive = is_positive(y[0])
    list_of_roots = []
    #for i in range(len(y)):
    for i in range(lenY):
        current_positive = is_positive(y[i])
        if current_positive != positive:
            list_of_roots.append(x[i - 1] + (x[i] - x[i - 1]) / (abs(y[i]) + abs(y[i - 1])) * abs(y[i - 1]))
            positive = not positive
    if len(list_of_roots) >= 2:
        if not return_as_dict:
            return abs(list_of_roots[-1] - list_of_roots[0])
        else:
            return {
                'fwhm': abs(list_of_roots[-1] - list_of_roots[0]),
                'x_range': list_of_roots,
            }
    else:
        raise Exception('Number of roots is less than 2!')

#****************************************************************************
def cal_coherence_length(x, y, shift=0.5, return_as_dict=False):
    """The function searches x-values (roots) where y=0 (after normalization to values between 0 and 1 and shifting the
    values down by 0.5 (default value)) based on linear interpolation, and calculates full width at half maximum (FWHM).

    :param x: an array of x values.
    :param y: an array of y values.
    :param shift: an optional shift to be used in the process of normalization (between 0 and 1).
    :param return_as_dict: if to return a dict with 'fwhm' and 'x_range'
    :return: a value of the FWHM or dictionary consisting of 'fwhm' and 'x_range'
    """

    def is_positive(num):
        return True if num > 0 else False

    #OC18112017 (making it work with standard Python lists / arrays)
    minY = min(y)
    maxY = max(y)
    maxY_ind = np.argmax(y)
    if(maxY == minY): raise Exception('coherence length can not be calculated')
    mult = 1./(maxY - minY)
    lenY = len(y)
    for i in range(lenY): y[i] = (y[i] - minY)*mult - shift

    i = maxY_ind
    for i in range(maxY_ind,lenY):
        if is_positive(y[i]) == False: 
            x_r = x[i - 1] + (x[i] - x[i - 1]) / (abs(y[i]) + abs(y[i - 1])) * abs(y[i - 1])
            break
        else:
            x_r = x[-1]
#            raise Exception('All y is larger than 0.5!')
    i = maxY_ind
    for i in reversed(range(maxY_ind + 1)):
        if is_positive(y[i]) == False: 
            x_l = x[i] + (x[i + 1] - x[i]) / (abs(y[i + 1]) + abs(y[i])) * abs(y[i])
            break
        else:
            x_l = x[0]
#           raise Exception('All y is larger than 0.5!')
    coh_len = x_r - x_l
    return coh_len

#****************************************************************************

def read_coherence_data_file(_fname):
    import uti_plot_com
    import uti_io
    import uti_math
    data, mode, allrange, arLabels, arUnits = uti_plot_com.file_load(_fname, 0, 0)
    #allrange, units = uti_plot_com.rescale_range(allrange, arUnits, 0, 0, 0)
    e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
    allrange = list(allrange)
    if nx == 1 and ny > 1:     
        allrange[3], allrange[4], allrange[5] = allrange[6], allrange[7], allrange[8]
    elif ny == 1 and nx > 1: 
        allrange[6], allrange[7], allrange[8] = allrange[3], allrange[4], allrange[5]
    allrange = tuple(allrange)
    e0, e1, ne, x0, x1, nx, y0, y1, ny = allrange
    data = np.array(data)

    xStep = (x1 - x0)/(nx - 1)
    yStep = (y1 - y0)/(ny - 1)
    inperpOrd = 1 #interpolation order to use (1 to 3)

    if _fname.split('.')[3] == '1': # horizontal coherent 
#        print('horizontal')
        y = np.linspace(y0, y1, ny)
        arCutY = array('d', [0]*ny)
        yy = y0
        for iy in range(ny):
            arCutY[iy] = uti_math.interp_2d(0, yy, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
            yy += yStep
        return y, arCutY

    if _fname.split('.')[3] == '2': # vertical coherent
#        print('vertical')
        x = np.linspace(x0, x1, nx)
        arCutX = array('d', [0]*nx)
        xx = x0
        for ix in range(nx):
            arCutX[ix] = uti_math.interp_2d(xx, 0, x0, xStep, nx, y0, yStep, ny, data, inperpOrd, 1, 0)
            xx += xStep
        return x, arCutX
#def Coherence_Length(_fname):
#    xd, yd = read_coherence_data_file(_fname)
#    coh_len = cal_coherence_length(xd, yd)

def Coherence_Length(_file_path_deg_coh1, _file_path_deg_coh2):
    if os.path.isfile(_file_path_deg_coh1):
        xd_1, yd_1 = read_coherence_data_file(_file_path_deg_coh1)
        coh_1 = cal_coherence_length(xd_1, yd_1)
    else:
        print('Warning! no Horizontal Coherent data file!')
        coh_1 = 1e99
    if os.path.isfile(_file_path_deg_coh2):
        xd_2, yd_2 = read_coherence_data_file(_file_path_deg_coh2)
        coh_2 = cal_coherence_length(xd_2, yd_2)
    else:
        print('Warning! no Vertical Coherent data file!')
        coh_2 = 1e99
    return coh_1, coh_2

def remove_file(_file_path):
    if os.path.isfile(_file_path):
        os.remove(_file_path)

