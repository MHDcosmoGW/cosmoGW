"""
hydro_bubbles.py is a Python routine that contains functions to study the 1D hydrodynamic
solutions of expanding bubbles produced from first-order phase transitions.

Author: Alberto Roper Pol
Created: 01/02/2023

Other contributors: Antonino Midiri, Simona Procacci

Main reference is:

Appendix A of RPPMC24 - A. Roper Pol, S. Procacci, A. S. Midiri, C. Caprini,
"Irrotational fluid perturbations from first-order phase transitions,"
in preparation

Other used references are:

EKNS10 - J. R. Espinosa, T. Konstandin, J. M. No, G. Servant,
"Energy Budget of Cosmological First-order Phase Transitions,"
JCAP 06 (2010) 028, arXiv:1004.4187

H16 - M. Hindmarsh, "Sound shell model for acoustic gravitational
wave production at a first-order phase transition in the early Universe,"
Phys.Rev.Lett. 120 (2018) 7, 071301, arXiv: 1608.04735

HH19 - M. Hindmarsh, M. Hijazi, "Gravitational waves from first order
cosmological phase transitions in the Sound Shell Model,"
JCAP 12 (2019) 062, arXiv:1909.10040

RPPC23 - A. Roper Pol, S. Procacci, C. Caprini,
"Characterization of the gravitational wave spectrum from sound waves within
the sound shell model," Phys. Rev. D 109, 063531 (2024), arXiv:2308.12943.
"""

import numpy as np
import matplotlib.pyplot as plt

### Reference values
cs2 = 1/3      # speed of sound
cols_ref = ['black', 'darkblue', 'blue', 'darkgreen', 'green',
            'purple', 'darkred', 'red', 'darkorange', 'orange', 'violet']

def Chapman_Jouget(alp):
    
    '''
    Chapman-Jouget is the wall velocity at which the relative speed behind the wall becomes
    that of the spped of sound. It is the limiting case separating detonations
    and supersonic deflagrations.
    
    Reference: Eq. 39 of EKNS10, Eq. B19 of HH19 has a typo
    
    Arguments:
        alp -- alpha at the + side of the wall (alpha_pl)
    '''

    return 1/np.sqrt(3)/(1 + alp)*(1 + np.sqrt(alp*(2 + 3*alp)))

def type_nucleation(vw, alp, cs2=cs2):
    
    '''
    Function that determines the type of bubble solution: subsonic deflagrations ('def'),
    supersonic deflagrations ('hyb'), or detonations ('det').
    
    Arguments:
        vw -- bubble wall speed
        alp -- phase transition parameter (based on enthalpy at the symmetric phase + side of the wall)
               Note that this alpha might be different than that at the temperature of nucleation
        cs2 -- square of the speed of sound (default 1/3)
        
    Returns:
        ty -- type of solution ('def', 'hyb', 'det')
    '''

    v_cj = Chapman_Jouget(alp)
    cs = np.sqrt(cs2)

    # check if vw is a list or a single value
    if not isinstance(vw, (list, tuple, np.ndarray)):

        if vw < cs: ty = 'def'
        elif vw < v_cj: ty = 'hyb'
        else: ty = 'det'

    else:

        ty = np.array(['hyb']*len(vw))
        ty[np.where(vw < cs)] = 'def'
        ty[np.where(vw > v_cj)] = 'det'

    return ty

#################### 1D HYDRO SOLUTIONS UNDER SPHERICAL SYMMETRY ####################

'''
The solutions are based in the 1D hydrodynamic descriptions of EKNS10 and HH19.
The details accompanying the code are provided in RPPC23 and RPPMC24 (appendix A)
'''

def Lor_mu(v, vw):
    
    '''
    Lorentz transform of the velocity v in the reference frame of the wall
    (with speed vw)
    
    mu(vw, v) = (vw - v)/(1 - vw*v)
    
    Reference: Eq. B12 of HH19
    '''

    return (vw - v)/(1 - v*vw)

### hydrodynamic shocks with no vacuum energy

def v_shock(xi):
    
    '''
    Function that computes the shock velocity at the - side of the
    shock. Computed from the v_+ v_- = 1/3 condition.
    '''
    
    vsh = (3*xi**2 - 1)/2/xi

    return vsh

def w_shock(xi):

    '''
    Function that computes the ratio of enthalpies w-/w+ across the
    shock. Computed from the v_+ v_- = 1/3 condition.
    '''

    wsh = (9*xi**2 - 1)/3/(1 - xi**2)

    return wsh

#### differential equation for the velocity radial profile

def xi_o_v(xi, v, cs2=cs2):
 
    '''
    Function that characterizes the 1d hydro equation under radial
    symmetry. It computes the value of dxi/dv used to solve the
    equation in compute_xi_from_v using RK4.

    Reference: Eq. 27 of EKNS10
        
    Arguments:
        xi -- self-similar r/t
        v -- 1d velocity profile
        cs2 -- square of the speed of sound (default 1/3)

    Returns:
        f -- value of f(xi, v) = dxi/dv
    '''

    gamma2 = 1/(1 - v**2)
    mu = Lor_mu(v, xi)
    f = xi*gamma2*(1 - xi*v)*(mu**2/cs2 - 1)/2/v

    return f

def compute_xi_from_v(v, xi0, cs2=cs2, shock=False):

    '''
    Function that computes the solution xi (v) using a 4th-order
    Runge-Kutta scheme. Since dv/dxi has a singularity, it is necessary
    to compute xi(v) and then invert each of the solutions to have the full
    dynamical solution. However, for the physical solution, computing v (xi)
    is more practical.

    Arguments:
        v -- velocity array
        xi0 -- position where boundary condition is known (it has to correspond
               to the first value in the velocity array)
        cs2 -- square of the speed of sound (default 1/3)
        shock -- option to stop the integration when a shock is found
                 (valid for deflagrations)
    '''

    xi = np.zeros(len(v)) + xi0
    sh = False
    indsh = -1

    for i in range(0, len(v) - 1):
        dv = v[i+1] - v[i]
        k1 = xi_o_v(xi[i], v[i])
        k2 = xi_o_v(xi[i] + dv*k1/2, .5*(v[i+1] + v[i]))
        k3 = xi_o_v(xi[i] + dv*k2/2, .5*(v[i+1] + v[i]))
        k4 = xi_o_v(xi[i] + dv*k3, v[i + 1])
        xi_new = xi[i] + 1/6*(k1 + 2*k2 + 2*k3 + k4)*dv
        if shock:
            xi_sh = xi_new
            v_sh = v_shock(xi_sh)
            if v[i + 1] < v_sh:
                xi[i + 1:] = xi_sh
                sh = True
                indsh = i
                break
        xi[i + 1] = xi_new
        if xi_new > 1: xi[i + 1] = 1

    return xi, sh, indsh

def compute_int_w(xi, v, cs2=cs2):
    
    '''
    Function that computes the integrand for the integration of dw/dxi
    (enthalpy) equation, as a function of the solution v(xi).
    
     Arguments:
        xi -- self-similar r/t
        v -- 1d velocity profile
        cs2 -- square of the speed of sound (default 1/3)
    '''

    return (1 + 1/cs2)/(1 - v**2)*Lor_mu(v, xi)

def compute_w(v, xi, cs2=cs2):

    '''
    Function that computes the enthalpy from the solution of the
    velocity profile.
    
     Arguments:
        xi -- self-similar r/t
        v -- 1d velocity profile
        cs2 -- square of the speed of sound (default 1/3)
    '''
    
    w = np.zeros(len(v)) + 1
    ss = 0
    for i in range(0, len(v) - 1):
        ff_ip = compute_int_w(xi[i + 1], v[i + 1], cs2=cs2)
        ff_i = compute_int_w(xi[i], v[i], cs2=cs2)
        ss += .5*(v[i + 1] - v[i])*(ff_ip + ff_i)
        w[i + 1] = np.exp(ss)
        
    return w

######### SOLVE FOR THE DIFFERENT TYPE OF BOUNDARY CONDITIONS ############

################# MATCHING CONDITIONS ACROSS DISCONTINUITIES #################

def vp_tilde_from_vm_tilde(vw, alpha, plus=True, sg='plus'):
    
    '''
    Function that computes the + (symmetric phase) or - (broken phase)
    velocity, defined in the wall reference frame, across the wall,
    as a function of the value of the velocity at the opposite side of
    the wall (via the matching conditions).
    
    Reference: Eqs. B6 and B7 of HH19
    
    Arguments:
        vw -- wall velocity, imposed at one of the sides of the wall
        alpha -- phase transition strength at the symmetric phase (can be different
                 than that at the nucleation temperature)
        plus -- option to consider positive (True) or negative (False) branch of the equation
                found from the matching conditions (default is True
        sg -- option to compute v+ from v- = vw if sg == 'plus' or to compute
              v- from v+ = vw if sg == 'minus' (default is 'plus')
              
    Returns:
        vp_vm -- v+ or v- from the value at the other side of the bubble wall
    '''
    
    if sg == 'plus':
        
        a1 = 1/3/vw + vw
        a2 = np.sqrt((1/3/vw - vw)**2 + 4*alpha**2 + 8/3*alpha)
        aa = .5/(1 + alpha)
        
    else:
        
        a1 = (1 + alpha)*vw + (1 - 3*alpha)/3/vw
        a2 = np.sqrt(((1 + alpha)*vw + (1 - 3*alpha)/3/vw)**2 - 4/3)
        aa = 1/2
        
    if plus: a = a1 + a2
    else: a = a1 - a2
    
    vp_vm = aa*a
    
    return vp_vm

def vplus_vminus(alpha, vw=1., ty='det', cs2=cs2):
    
    '''
    Function that returns \tilde v_+ and \tilde v_- for the different
    type of solutions (deflagrations, detonations and hybrids).
    This allows to give the boundary conditions corresponding to each
    of the solutions.
    
    Arguments:
        alpha -- value of alpha at the + side of the wall
        vw -- wall velocity
        cs2 -- speed of sound squared (default is 1/3)
        ty -- type of solution
        
    Returns:
        vplus, vminus -- + and - velocities expressed in the wall
                            reference of frame
    '''

    cs = np.sqrt(cs2)
    
    if not isinstance(ty, (list, tuple, np.ndarray)):

        if ty == 'det':
            vplus = vw
            vminus = vp_tilde_from_vm_tilde(vw, alpha, plus=True, sg='minus')
        if ty == 'def':
            vminus = vw
            vplus = vp_tilde_from_vm_tilde(vw, alpha, plus=False, sg='plus')
        if ty == 'hyb':
            vminus = cs
            vplus = vp_tilde_from_vm_tilde(cs, alpha, plus=False, sg='plus')
   
    else:

        vplus = np.zeros(len(vw))
        vminus = np.zeros(len(vw))
        inds_det = np.where(ty == 'det')
        inds_def = np.where(ty == 'def')
        inds_hyb = np.where(ty == 'hyb')
        vplus[inds_det] = vw[inds_det]
        vminus[inds_det] = vp_tilde_from_vm_tilde(vw, alpha, plus=True, sg='minus')[inds_det]
        vplus[inds_def] = vp_tilde_from_vm_tilde(vw, alpha, plus=False, sg='plus')[inds_def]
        vminus[inds_def] = vw[inds_def]
        vplus[inds_hyb] = vp_tilde_from_vm_tilde(cs, alpha, plus=False, sg='plus')
        vminus[inds_hyb] = cs

    return vplus, vminus

################## function that computes the detonation part of the solutions ##################
def det_sol(v0, xi0, cs2=cs2, Nxi=1000, zero_v=-4):
    
    '''
    Function that computes a detonation solution with boundary condition
    v0 at xi0
    
    Arguments:
        v0 -- value of v0 at the boundary
        xi0 -- position of the boundary
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi
        zero_v -- reference zero velocity (default 1e-4)
        
    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
    '''
    
    # compute solution from initial condition v = v0 at xi = vw
    # until v reduces to 4 orders of magntiude below value v minus
    cs = np.sqrt(cs2)

    # included option to initialize with multiple vws
    if not isinstance(xi0, (list, tuple, np.ndarray)):

        vs = np.logspace(np.log10(v0), np.log10(v0) + zero_v, Nxi)
        xis, sh, indsh = compute_xi_from_v(vs, xi0, cs2=cs2, shock=False)
        xi_sh = xis[indsh]
        ws = compute_w(vs, xis, cs2=cs2)

        inds_sort = np.argsort(xis)
        xis = xis[inds_sort]
        vs = vs[inds_sort]
        ws = ws[inds_sort]
        
    else:

        xis = np.zeros((len(xi0), Nxi))
        vs = np.zeros((len(xi0), Nxi))
        ws = np.zeros((len(xi0), Nxi))

        for i in range(0, len(xi0)):

            vs[i, :] = np.logspace(np.log10(v0), np.log10(v0) + zero_v, Nxi)
            xis[i, :], _, _ = compute_xi_from_v(vs[i, :], xi0[i], cs2=cs2, shock=False)
            ws[i, :] = compute_w(vs[i, :], xis[i, :], cs2=cs2)
            
            inds_sort = np.argsort(xis[i, :])
            xis[i, :] = xis[i, inds_sort]
            vs[i, :] = vs[i, inds_sort]
            ws[i, :] = ws[i, inds_sort]

    return xis, vs, ws

######## function that computes the deflagration part of the solutions ###########
def def_sol(v0, xi0, cs2=cs2, Nxi=1000, shock=True, zero_v=-4):
    
    '''
    Function that computes a deflagration solution with boundary condition
    v0 at xi0
    
    Arguments:
        v0 -- value of v0 at the boundary
        xi0 -- position of the boundary
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi
        shock -- possibility to stop the calculation once a shock is formed
                 (default is True)
        zero_v -- reference zero velocity (default 1e-4)

    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
        xi_sh, sh -- position of the shock and boolean which becomes
                        True if a shock forms
    '''

    vs = np.logspace(np.log10(v0), np.log10(v0) + zero_v, Nxi)
    cs = np.sqrt(cs2)
    xi_sh = cs
    sh = False
    v_sh = 0
    
    if shock:
        xiss, sh, indsh = compute_xi_from_v(vs, xi0, cs2=cs2, shock=True)
        xi_sh = xiss[indsh]
        if not sh: xi_sh = cs
        v_sh = v_shock(xi_sh)
        xis = np.linspace(xi0, xi_sh, Nxi + 1)
        xis = xis[:Nxi-1]
        vs = np.interp(xis, xiss, vs)
        vs = np.append(vs, v_sh)
        xis = np.append(xis, xi_sh)

    else:
        xis, sh, indsh = compute_xi_from_v(vs, xi0, cs2=cs2, shock=False)
        xi_sh = xis[indsh]
    
    ws = compute_w(vs, xis, cs2=cs2)    
    w_sh = w_shock(xi_sh)
    ws = ws*w_sh/ws[-1]
    
    return xis, vs, ws, xi_sh, sh

def compute_def(vw=0.5, alpha=0.263, cs2=cs2, Nxi=1000, shock=True):
    
    '''
    Function that computes the solutions for a subsonic deflagration
    1d profile given vw and alpha, using def_sol
    
    Arguments:
        vw -- wall velocity (default is 0.5)
        alpha -- strength of the phase transition (default is 0.263)
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi
        shock -- possibility to stop the calculation once a shock is formed
                 (default is True)

    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
        xi_sh, sh -- position of the shock and boolean which becomes
                        True if a shock forms
        w_pl, w_m -- plus and minus values of the enthalpies
                        across the bubble
    '''
    
    ## relative velocity at + is computed from \tilde v- = \xi_w
    vrels, _ = vplus_vminus(alpha, vw=vw, ty='def')
    # Lorentz boosted v plus
    vpl = Lor_mu(vrels, vw)
    
    xis, vs, ws, xi_sh, sh = def_sol(vpl, vw, cs2=cs2, Nxi=Nxi, shock=shock)
    
    # values at both sides of the bubble wall
    w_pl = ws[0]
    w_m = w_pl*vrels/(1 - vrels**2)/vw*(1 - vw**2)

    return xis, vs, ws, xi_sh, sh, w_pl, w_m

def compute_hyb(vw=0.7, alpha=0.052, cs2=cs2, Nxi=1000, shock=True):
    
    '''
    Function that computes the solutions for a supersonic deflagration
    1d profile given vw and alpha, using det_sol and def_sol
    
    Arguments:
        vw -- wall velocity (default is 0.7)
        alpha -- strength of the phase transition (default is 0.052)
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi
        shock -- possibility to stop the calculation once a shock is formed
                 (default is True)

    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
        xi_sh, sh -- position of the shock and boolean which becomes
                        True if a shock forms
        w_pl, w_m -- plus and minus values of the enthalpies
    '''
    
    cs = np.sqrt(cs2)
    
    ## relative velocity at + is computed from \tilde v- = cs
    vrels, _ = vplus_vminus(alpha, cs2=cs2, ty='hyb')
    vpl = Lor_mu(vrels, vw)
    vm = Lor_mu(cs, vw)
    
    # compute deflagration solution
    xis, vs, ws, xi_sh, sh = def_sol(vpl, vw, cs2=cs2, Nxi=int(Nxi/2), shock=shock)
    
    # compute detonation solution
    xis2, vs2, ws2 = det_sol(vm, vw, cs2=cs2, Nxi=int(Nxi/2))
    # ratio of w+ over w- across the bubble wall
    w_pl = ws[0]
    w_m = w_pl*vrels/(1 - vrels**2)*(1 - cs2)/cs
    ws2 *= w_m
    
    xis = np.append(xis2, xis)
    vs = np.append(vs2, vs)
    ws = np.append(ws2, ws)

    return xis, vs, ws, xi_sh, sh, w_pl, w_m

def compute_det(vw=0.77, alpha=0.091, cs2=cs2, Nxi=1000):

    '''
    Function that computes the solutions for a detonation 1d profile
    given vw and alpha, using det_sol
    
    Arguments:
        vw -- wall velocity (default is 0.77)
        alpha -- strength of the phase transition (default is 0.091)
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi

    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
        w_pl, w_m -- plus and minus values of the enthalpies
                        across the bubble
    '''
    
    ## relative velocity at - is computed from \tilde v+ = \xi_w
    _, vrels = vplus_vminus(alpha, vw=vw, ty='det')
    # Lorentz boosted v minus
    vm = Lor_mu(vrels, vw)
    w_m = vw/(1 - vw**2)/vrels*(1 - vrels**2)
    w_pl = 1
    
    xis, vs, ws = det_sol(vm, vw, cs2=cs2, Nxi=Nxi)
    ws *= w_m
    
    # no shock is formed in detonations, so xi_sh is set to vw
    # and sh to False
    xi_sh = vw
    sh = False

    return xis, vs, ws, xi_sh, sh, w_pl, w_m

def compute_alphan(vw=0.5, alpha_obj=0.263, tol=1e-4, cs2=cs2, quiet=False,
                   max_it=30, meth=1, Nxi=1000, ty='def'):
    
    '''
    Function that computes the value of \alpha_+ corresponding to alpha.
    It requires to compute the 1d profile of w and then iteratively look for
    the value of alpha_+ that gives the correct alpha.
    
    Arguments:
        vw -- wall velocity
        alpha_obj -- value of alpha defined at the nucleation temperature that
                        wants to be reached
        tol -- relative tolerance to consider convergence (default is 1e-4)
        cs2 -- speed of sound squared (default is 1/3)
        quiet -- option to avoid printing some debug options (default is False)
        max_it -- maximum number of allowed iterations (default is 30)
        meth -- method on the Newton-Raphson update (2 options)
        Nxi -- number of points in xi discretization (default is 1000)
        ty -- type of solution (options are def or hyb)
        
    Returns:
        xis0, vvs0, wws0 -- arrays of xi, velocity and enthalpy of the converged
                               solutions
        xi_sh, sh, w_pl, w_m -- shock position, boolean if shock has formed,
                                    plus and minus enthalpies
        alpha_n -- converged (or not) alpha at nucleation
        alp_plus -- value of alpha_+ leading to alpha_obj
        conv -- boolean determining if the algorithm has converged
    '''
    
    alps_tst = np.logspace(-2, 2, 10000)
    vCJs_tst = Chapman_Jouget(alps_tst)
    # first guess
    alp_plus = alpha_obj

    j = 0
    conv = False

    while not conv and j < max_it:
        
        j += 1

        if ty == 'hyb':
            xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m = \
                compute_hyb(vw=vw, alpha=alp_plus, cs2=cs2, Nxi=Nxi, shock=True)

        if ty == 'def':
            xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m = \
                compute_def(vw=vw, alpha=alp_plus, cs2=cs2, Nxi=Nxi, shock=True)
                
        alpha_n = alp_plus*w_pl
        
        if abs(alpha_n - alpha_obj)/alpha_obj < tol: conv = True

        else:

            if meth==1: alp_plus = alpha_obj/w_pl
            if meth==2: alp_plus += (alpha_obj - alpha_n)/w_pl
        if not quiet:
            print('iteration', j, 'alpha', alpha_n)
            print('iteration', j, 'new guess', alp_plus)

    print(j, 'iterations for vw=', vw,' and alpha= ',alpha_obj)
    print('alpha:', alpha_n, ', alpha_+:', alp_plus)
        
    return xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m, alpha_n, alp_plus, conv

def compute_profiles_vws(alpha, vws=[], cs2=cs2, Nvws=20, Nxi=10000, Nxi2=10000, plot=True, plot_v='v', cols=[],
                         alphan=True, quiet=True, tol=1e-5, max_it=30, ls='solid', alp=1., lam=False, meth=1,
                         legs=False, fs_lg=14, st_lg=2, eff=False, save=False, dec_vw=1, ress='results/1d_profiles',
                         strs_vws=[], str_alp=[]):

    '''
    Function that computes the velocity and enthalpy profiles for a given
    alpha (at nucleation T, not alpha+) and a range of wall velocities.
    
    Arguments:
        alpha -- nucleation T alpha
        vws -- range of wall velocities
        cs2 -- square of the speed of sound (default is 1/3)
        Nvws -- number of wall velocities (if vws is not given)
        Nxi -- number of discretization points in xi where the profiles are
               computed
        Nxi2 -- number of discretization points in xi from the constructed
                profiles from 0 to 1
        plot -- option to plot the resulting 1d profiles
        plot_v -- choice of plotting ('v' for velocity, 'w' for enthalpy, 'both' for both)
        cols -- array with colors (cols_ref by default)
        alphan -- option to identify if the input alpha is at the nucleation temperature
                  (if True) or alpha+ (if False), default is True
        quiet -- option to avoid printing debugging information (default is True)
        max_it -- maximum number of iterations to find alpha_+
        ls -- line styles
        alp -- opacity of the plots
        lam -- option to compute energy perturbations lambda instead of enthalpy (default is False,
               so enthalpy)
        tol -- tolerance of the relative error to consider convergence of alpha_+ has been reached
        legs -- legends to be included (default is False)
        save -- option to save in a file the results of the 1d profiles (default is False)
        ress -- directory where to save the files (default is 'results/1d_profiles')
        eff -- option to compute the efficiency factors kappa and omega (default is False)
    '''

    if save: import pandas as pd

    if vws == []: vws = np.linspace(0.1, .99, Nvws)
    vCJ = Chapman_Jouget(alp=alpha)
    cs = np.sqrt(cs2)
    xis = np.linspace(0, 1, Nxi2)
    vvs = np.zeros((len(vws), len(xis)))
    wws = np.zeros((len(vws), len(xis))) + 1
    alphas_n = np.zeros(len(vws))
    conv = np.zeros(len(vws)) + 1
    kappas = np.zeros(len(vws))
    omegas = np.zeros(len(vws))
    shocks = np.zeros(len(vws))
    xi_shocks = np.zeros(len(vws))
    wms = np.zeros(len(vws))
    
    if plot_v == 'both' and plot:
        plt.figure(1)
        plt.figure(2)
    
    if cols == []: cols = cols_ref

    for i in range(0, len(vws)):
        
        # determine type of solution
        ty = type_nucleation(vws[i], alpha, cs2=cs2)
        
        j = i%11
        if ty == 'def':
            
            ## iteratively compute the real alpha_+ leading to alpha
            if alphan:
                xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m, alpha_n, alp_plus, conv[i] = \
                        compute_alphan(vw=vws[i], alpha_obj=alpha, tol=tol, cs2=cs2, meth=meth,
                                       quiet=quiet, max_it=max_it, Nxi=Nxi, ty='def')
                alphas_n[i] = alp_plus

            ## otherwise, alpha given is assumed to be alpha_+
            else:
                xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m = \
                        compute_def(vw=vws[i], alpha=alpha, cs2=cs2, Nxi=Nxi, shock=True)
                alphas_n[i] = alpha*w_pl
            
            inds = np.where((xis >= vws[i])*(xis <= xi_sh))[0]
            inds2 = np.where(xis < vws[i])[0]
            wws[i, inds2] = w_m            

        elif ty == 'hyb':

            ## iteratively compute the real alpha_+ leading to alpha
            if alphan:
                xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m, alpha_n, alp_plus, conv[i] = \
                        compute_alphan(vw=vws[i], alpha_obj=alpha, tol=tol, cs2=cs2, ty='hyb',
                                       meth=meth, quiet=quiet, max_it=max_it)
                alphas_n[i] = alp_plus

            ## otherwise, alpha given is assumed to be alpha_+
            else:
                xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m = \
                        compute_hyb(vw=vws[i], alpha=alpha, cs2=cs2,
                                    Nxi=Nxi, shock=True)
                alphas_n[i] = alpha*w_pl

            inds = np.where((xis >= cs)*(xis <= xi_sh))[0]
            inds2 = np.where(xis < cs)[0]
            wws[i, inds2] = wws0[0]

        else:
            
            xis0, vvs0, wws0, xi_sh, sh, w_pl, w_m = \
                    compute_det(vw=vws[i], alpha=alpha, cs2=cs2,
                                Nxi=Nxi)

            inds = np.where((xis >= cs)*(xis <= vws[i]))[0]
            inds2 = np.where(xis < cs)[0]
            wws[i, inds2] = wws0[0]
            alphas_n[i] = alpha
        
        vvs[i, inds] = np.interp(xis[inds], xis0, vvs0)
        wws[i, inds] = np.interp(xis[inds], xis0, wws0)
        shocks[i] = sh
        xi_shocks[i] = xi_sh
        wms[i] = w_m

        # compute efficiency of energy density production
        if eff:
            kappas[i], omegas[i] = kappas_from_prof(vws[i], alpha,
                                                    xis, wws[i, :], vvs[i, :])
      
        # compute mean energy density from enthalpy if lam is True
        if lam:
            if alphan: alp_lam = alpha
            else: alp_lam = alphas_n[i]
            wws[i, :] = w_to_lam(xis, wws[i, :], vws[i], alp_lam)

        if plot:

            if st_lg == 1: str_lg = r'$\xi_w=%.1f$'%vws[i]
            if st_lg == 2: str_lg = r'$\xi_w=%.2f$'%vws[i]
            if plot_v=='v': plt.plot(xis, vvs[i, :], color=cols[j], ls=ls, alpha=alp,
                                     label=str_lg)
            if plot_v=='w': plt.plot(xis, wws[i, :], color=cols[j], ls=ls, alpha=alp,
                                     label=str_lg)
            if plot_v=='both':
                plt.figure(1)
                plt.plot(xis, vvs[i, :], color=cols[j], ls=ls, alpha=alp, label=str_lg)
                plt.figure(2)
                plt.plot(xis, wws[i, :], color=cols[j], ls=ls, alpha=alp, label=str_lg)
                
        if save and alphan:

            # option to save the results on files (one file for each alpha and wall velocity)
            # save is only possible when the input alpha is that at the nucleation temperature (alphan = True)
            df = pd.DataFrame({'alpha': alpha*xis**0, 'xi_w': vws[i]*xis**0, 'xi': xis, 'v': vvs[i, :],
                               'w': wws[i, :], 'alpha_pl': alphas_n[i]*xis**0, 'shock': shocks[i]*xis**0,
                               'xi_sh': xi_shocks[i]*xis**0, 'wm': wms[i]*xis**0})
            # save file
            if str_alp == []:
                str_alp = '%s'%(alpha)
                str_alp = '0' + str_alp
                str_alp = str_alp[2:]
            if strs_vws == []:
                str_vws = '%s'%(np.round(vws[i], decimals=dec_vw))
                str_vws = str_vws[2:]
            else: str_vws = strs_vws[i]
            file_dir = ress + '/alpha_%s_vw_%s.csv'%(str_alp, str_vws)
            df.to_csv(file_dir)
            print('results of 1d profile saved in ', file_dir)

    if plot:

        if plot_v == 'v' or plot_v == 'both':
            if plot_v == 'both': plt.figure(1)
            plt.ylim(-.05, 1.05)
            plt.ylabel(r'$ v_{\rm ip} (\xi)$')
        if plot_v == 'w' or plot_v == 'both':
            if plot_v == 'both': plt.figure(2)
            plt.ylim(0, 5)
            if lam: plt.ylabel(r'$ \lambda_{\rm ip} (\xi)$')
            else: plt.ylabel(r'$ w(\xi)$')
        l = [1]
        if plot_v == 'both': l = [1, 2]
        for j in l: 
            plt.figure(j)
            plt.xlim(0, 1)
            plt.vlines(cs, -5, 30, color='black', ls='dashed', lw=1)
            plt.vlines(vCJ, -5, 30, color='black', ls='dashed', lw=1)
            plt.xlabel(r'$\xi$')
            if legs: plt.legend(fontsize=fs_lg)
        
    if eff:
        return xis, vvs, wws, alphas_n, conv, shocks, xi_shocks, wms, kappas, omegas
    else:
        return xis, vvs, wws, alphas_n, conv, shocks, xi_shocks, wms

################### COMPUTING EFFICIENCIES FROM 1D PROFILES ###################

# def kappas_from_prof(vw, alpha, xis, ws, vs, more=False, cs2=cs2, sw=False):

#     '''
#     Function that computes the kinetic energy density efficiency kappa
#     and thermal factor omega from the 1d profiles.
#     '''

#     kappa = 4/vw**3/alpha*np.trapz(xis**2*ws/(1 - vs**2)*vs**2, xis)
#     omega = 3/vw**3/alpha*np.trapz(xis**2*(ws - 1), xis)

#     return kappa, omega

# adapted from GW_fopt, combined kappa and kappas, to be tested
def kappas_Esp(vw, alp, cs2=cs2):

    """"
    Function that computes the efficiency in converting vacuum to
    kinetic energy density for detonations, deflagrations and hybrids.

    Uses the semiempirical fits from EKNS10, appendix A.

    Arguments:
        vw -- value or array of wall velocities
        alp -- strength of the phase transition at the nucleation temperature
        cs2 -- square of the speed of sound (default is 1/3)

    Returns:
        kappa_def, kappa_det, kappa_hyb: efficiency of kinetic energy production
            as a fraction of alpha/(1 + alpha) assuming the profile is a deflagration
            ('def'), a detonation ('det') and a hybrid supersonic deflagration ('hyb')
    """

    cs = np.sqrt(cs2)
    ty = type_nucleation(vw, alp, cs2=cs2)
    multi = False

    if isinstance(alp, (list, tuple, np.ndarray)):
        alp, vw = np.meshgrid(alp, vw, indexing='ij')
        multi = True

    # kappa at vw << cs
    kapA = vw**(6/5)*6.9*alp/(1.36 - 0.037*np.sqrt(alp) + alp)
    # kappa at vw = cs
    kapB = alp**(2/5)/(0.017 + (0.997 + alp)**(2/5))
    # kappa at vw = cJ (Chapman-Jouget)
    kapC = np.sqrt(alp)/(0.135 + np.sqrt(0.98 + alp))
    v_cj = Chapman_Jouget(alp)

    # kappa at vw -> 1    
    kapD = alp/(0.73 + 0.083*np.sqrt(alp) + alp)
    kappa_def = cs**(11/5)*kapA*kapB/((cs**(11/5) - vw**(11/5))*kapB + vw*cs**(6/5)*kapA)
    kappa_det = (v_cj - 1)**3*(v_cj/vw)**(5/2)*kapC*kapD
    den_kappa = ((v_cj - 1)**3 - (vw - 1)**3)*v_cj**(5/2)*kapC + (vw - 1)**3*kapD
    kappa_det = kappa_det/den_kappa

    ddk = -.9*np.log10(np.sqrt(alp)/(1 + np.sqrt(alp)))
    kappa_hyb = kapB + (vw - cs)*ddk + ((vw - cs)/(v_cj - cs))**3*(kapC - kapB - (v_cj - cs)*ddk)

    if isinstance(vw, (list, tuple, np.ndarray)):

        inds_def = np.where(vw <= cs)[0]
        inds_hyb = np.where((vw > cs)*(vw <= v_cj))[0]
        inds_det = np.where(vw > v_cj)[0]

        kap = np.zeros(np.shape(alpha))
        # subsonic deflagrations
        kap[inds_def] = kappa_def[inds_def]
        # supersonic detonations
        kap[inds_det] = kappa_det[inds_det]
        # supersonic deflagrations (hybrid)
        kap[inds_hyb] = kappa_hyb[inds_hyb]

    else:

        # subsonic deflagrations
        if vw < cs: kap =  kappa_def
        # supersonic detonations
        elif vw > v_cj: kap =  kappa_det
        # supersonic deflagrations (hybrid)
        else: kap = kappa_hyb

    return kap

######################### COMPUTING DIAGNOSTIC PROFILES #########################
    
def w_to_lam(xis, ws, vw, alpha_n):
    
    inds = np.where(xis < vw)[0]
    lam = 3/4*(ws - 1)
    lam[inds] -= 3/4*alpha_n
    
    return lam

########## COMPUTING FUNCTIONS RELEVANT FOR VELOCITY SPECTRAL DENSITY ##########

'''
Main reference is: RPPMC24
Other used references are: H16, HH19, RPPC23
'''

#### f' and l functions
    
def fp_z(xi, vs, z, ls=[], multi=True, quiet=False, lz=False, gpz=False):
    
    '''
    Function that computes the functions f'(z) and l(z) that appears in the Fourier
    transform of the velocity and enthalpy perturbations fields of an expanding bubble.
    '''

    xi_ij, z_ij = np.meshgrid(xi[1:], z, indexing='ij')
    zxi_ij = z_ij*xi_ij
    j1_z = np.sin(zxi_ij)/zxi_ij**2 - np.cos(zxi_ij)/zxi_ij
    # avoid division by zero
    j1_z[np.where(zxi_ij == 0)] = 0

    if lz:
        if ls == []:
            print('if lz is chosen you need to provide a l(xi) profile')
            lz = False
        j0_z = np.sin(zxi_ij)/zxi_ij
        j0_z[np.where(zxi_ij == 0)] = 1

    if multi:
        Nvws = np.shape(vs)[0]
        fpzs = np.zeros((Nvws, len(z)))
        if lz: lzs = np.zeros((Nvws, len(z)))

        for i in range(0, Nvws):
            v_ij, z_ij = np.meshgrid(vs[i, 1:], z, indexing='ij')
            fpzs[i, :] = -4*np.pi*np.trapz(j1_z*xi_ij**2*v_ij, xi[1:], axis=0)
            if lz:
                l_ij, z_ij = np.meshgrid(ls[i, 1:], z, indexing='ij')
                lzs[i, :] = 4*np.pi*np.trapz(j0_z*xi_ij**2*l_ij, xi[1:], axis=0)
                
            if not quiet: print('vw ', i + 1, '/', Nvws, ' computed')

    else:

        v_ij, z_ij = np.meshgrid(vs[1:], z, indexing='ij')
        fpzs = -4*np.pi*np.trapz(j1_z*xi_ij**2*v_ij, xi[1:], axis=0)
        if lz:
            l_ij, z_ij = np.meshgrid(ls[1:], z, indexing='ij')
            lzs = 4*np.pi*np.trapz(j0_z*xi_ij**2*l_ij, xi[1:], axis=0)

    if lz:
        return fpzs, lzs
    else:
        return fpzs

############################### NOT PUBLIC ON GITHUB ##############################

def compute_vws(vws, alphapl, shocks, xifs, ty, cs2=cs2):

    '''
    Function that computes xif and xib for a range of wall velocities
    and alpha+
    '''

    vpls = np.zeros(len(vws))
    vms = np.zeros(len(vws))
    vfs = np.zeros(len(vws))
    wfs = np.zeros(len(vws))

    cs = np.sqrt(cs2)

    for i in range(0, len(vws)):

        vplus, vminus = vplus_vminus(alphapl[i], vw=vws[i], cs2=cs2, ty=ty[i])
        vpls[i] = Lor_mu(vplus, vws[i])
        vms[i] = Lor_mu(vminus, vws[i])

        if shocks[i]:
            vfs[i] = v_shock(xifs[i])
            wfs[i] = w_shock(xifs[i])

    return vpls, vms, vfs, wfs

################### COMPUTING EFFICIENCIES FROM 1D PROFILES ###################

def kappas_from_prof(vw, alpha, xis, ws, vs, more=False, cs2=cs2, sw=False):

    '''
    Function that computes the kinetic energy density efficiency kappa
    and thermal factor omega from the 1d profiles.
    '''

    kappa = 4/vw**3/alpha*np.trapz(xis**2*ws/(1 - vs**2)*vs**2, xis)
    omega = 3/vw**3/alpha*np.trapz(xis**2*(ws - 1), xis)

    if more:

        kappasv = 4/vw**3/alpha*np.trapz(xis**2*vs**2, xis)
        kappas_sub = 4/vw**3/alpha*np.trapz(xis**2*ws*vs**2, xis)

        ## this part to be made public after paper on sound-wave efficiency is ready
        if sw:
            lams = w_to_lam(xis, ws, vw, alpha)
            kappas_sw = 2/vw**3/alpha*np.trapz(xis**2*(vs**2 + cs2*lams[i, :]**2), xis)

    return kappa, omega

####
#### function that computes the 1d solutions from an initial condition v0, xi0 ####
#### removed from the public version, attempt to unify the functions def_sol and det_sol
#### work in progress, one needs to check whether this is more efficient
####

def oned_sol(v0, xi0, cs2=cs2, Nxi=1000, zero_v=-4, ty='det', shock=True):
    
    '''
    Function that computes a detonation or a deflagration
    solution with boundary condition v0 at xi0
    
    Arguments:
        v0 -- value of v0 at the boundary
        xi0 -- position of the boundary
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi
        zero_v -- reference zero velocity (default 1e-4)
        
    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
    '''

    # included option to initialize with multiple vws and types
    if not isinstance(xi0, (list, tuple, np.ndarray)):

        if ty == 'det': shock = False
        vs = np.logspace(np.log10(v0), np.log10(v0) + zero_v, Nxi)
        xis, sh, indsh = hb.compute_xi_from_v(vs, xi0, cs2=cs2, shock=shock)
        ws = hb.compute_w(vs, xis, cs2=cs2)

        # correct enthalpy for deflagrations using the value of w_sh
        if sh:
            xi_sh = xis[indsh]
            w_sh = hb.w_shock(xi_sh)
            ws = ws*w_sh/ws[indsh]
        else: xi_sh = max(xi0, np.sqrt(cs2))

        # order arrays in ascending xi
        inds_sort = np.argsort(xis)
        xis = xis[inds_sort]
        vs = vs[inds_sort]
        ws = ws[inds_sort]
        
        # cut solution above xi_sh if shock is formed
        if sh:
            xiss = np.linspace(xi0, xi_sh, Nxi)
            vs = np.interp(xiss, xis, vs)
            ws = np.interp(xiss, xis, ws)
            xis = xiss

    else:

        xis = np.zeros((len(xi0), Nxi))
        vs = np.zeros((len(xi0), Nxi))
        ws = np.zeros((len(xi0), Nxi))
        xi_sh = np.zeros(len(xi0))
        sh = np.zeros(len(xi0))

        for i in range(0, len(xi0)):

            if ty[i] == 'det': shock = False
            vs[i, :] = np.logspace(np.log10(v0[i]), np.log10(v0[i]) + zero_v, Nxi)
            xis[i, :], sh[i], indsh = hb.compute_xi_from_v(vs[i, :], xi0[i], cs2=cs2, shock=shock)
            ws[i, :] = hb.compute_w(vs[i, :], xis[i, :], cs2=cs2)

            if sh[i]:
                xi_sh[i] = xis[i, indsh]
                w_sh = hb.w_shock(xi_sh[i])
                ws[i, :] = ws[i, :]*w_sh/ws[i, indsh]
            else: xi_sh[i] = max(xi0[i], np.sqrt(cs2))

            # order arrays in ascending xi
            inds_sort = np.argsort(xis[i, :])
            xis[i, :] = xis[i, inds_sort]
            vs[i, :] = vs[i, inds_sort]
            ws[i, :] = ws[i, inds_sort]

            # cut solution above xi_sh if shock is formed
            if sh[i]:
                xiss = np.linspace(xi0[i], xi_sh[i], Nxi)
                vs[i, :] = np.interp(xiss, xis[i, :], vs[i, :])
                ws[i, :] = np.interp(xiss, xis[i, :], ws[i, :])
                xis[i, :] = xiss

    return xis, vs, ws, xi_sh, sh

####
#### function that attempts to combine compute_def, compute_hyb, compute_det in one
#### single function to reduce code (at the moment, it is not more efficient in time)
#### and not as robust, so it is still work in progress
####

def compute_profile(vw, alpha, cs2=cs2, Nxi=1000, shock=True):
    
    '''
    Function that computes the 1d profiles for a given vw and alpha
    values.
    
    Arguments:
        vw -- wall velocity
        alpha -- strength of the phase transition
        cs2 -- speed of sound squared (default is 1/3)
        Nxi -- number of discretization points in xi
        shock -- possibility to stop the calculation once a shock is formed
                 (default is True)

    Returns:
        xis -- array of xi
        vs -- array of 1d velocities
        ws -- array of 1d enthalpies
        xi_sh, sh -- position of the shock and boolean which becomes
                     True if a shock forms (if a shock does not form xi_sh is
                     taken to be cs for subsonic deflagrations and vw for detonations)
        w_pl, w_m -- plus and minus values of the enthalpies
                        across the bubble
        w0 -- value of the enthalpy at the interior of the bubble (it is equal to w_m
                     for deflagrations and detonations but can be different for hybrids)
    '''
    
    ty = type_nucleation(vw, alpha, cs2=cs2)
    cs = np.sqrt(cs2)
    
    if ty == 'hyb': Nxi = int(Nxi/2)
    
    if ty == 'def' or ty == 'hyb':
        
        ## relative velocity at + is computed from \tilde v- = \xi_w
        vrels, _ = vplus_vminus(alpha, vw=vw, cs2=cs2, ty=ty)
        # Lorentz boosted v plus
        vpl = Lor_mu(vrels, vw)
        if ty == 'hyb': vm = Lor_mu(cs, vw)
        
        # compute deflagration solution
        # xis, vs, ws, xi_sh, sh = def_sol(vpl, vw, cs2=cs2, Nxi=Nxi, shock=shock)
        
        # values at both sides of the bubble wall
        w_pl = ws[0]
        w_m = w_pl*vrels/(1 - vrels**2)/vw*(1 - vw**2)
    
    if ty == 'det':
        
        ## relative velocity at - is computed from \tilde v+ = \xi_w
        _, vrels = vplus_vminus(alpha, vw=vw, ty='det')
        # Lorentz boosted v minus
        vm = Lor_mu(vrels, vw)
        w_m = vw/(1 - vw**2)/vrels*(1 - vrels**2)
        w_pl = 1

        xis, vs, ws = det_sol(vm, vw, cs2=cs2, Nxi=Nxi)
        ws *= w_m

        # no shock is formed in detonations, so xi_sh is set to vw
        # and sh to False
        xi_sh = vw
        sh = False

    # compute detonation solution
    xis2, vs2, ws2 = det_sol(vm, vw, cs2=cs2, Nxi=int(Nxi/2))
    # ratio of w+ over w- across the bubble wall
    w_pl = ws[0]
    w_m = w_pl*vrels/(1 - vrels**2)*(1 - cs2)/cs
    ws2 *= w_m
    
    xis = np.append(xis2, xis)
    vs = np.append(vs2, vs)
    ws = np.append(ws2, ws)

def x_prof(xis, ws, vs, vw):
    
    '''
    Function to compute the 1d profile of the variable Xi = \sqrt(w) gamma vi
    '''
    
    xpr = vs*np.sqrt(ws/(1 - vs**2))
    
    return xpr

################# ANALYTICAL CALCULATIONS FOR THE FUNCTION f'(z) #################

'''
Main reference is: RPPMC24
'''

def fp_z0(xi, v):
    
    '''
    Function that computes the limit z -> 0 of the function f'(z)/z
    '''

    fpz0 = -4*np.pi/3*np.trapz(v*xi**3, xi)

    return fpz0

##############      f'(z) function for toy models      ##############

######################### Define toy models #########################

def toy_profile(xi, v=1, xif=1, xib=0, ty='const'):

    '''
    Function that returns the velocity profile as a constant or linear
    toy profile if the solution is a deflagration, hybrid or detonation.
    '''

    if ty == 'const': vtoy = xi**0
    elif ty == 'lin_inc': vtoy = (xi - xib)/(xif - xib)
    elif ty == 'lin_dec': vtoy =  (-xi + xif)/(xif - xib)
    else:
        print('Options for toy profile are const, lin_inc, and lin_dec')
        vtoy = 0.

    return vtoy*v

def fp_zinf(z, vw, Dvw, sh=False, vsh=1., xsh=1.):
    
    '''
    Function that computes the limit z -> infinity of the function f'(z)/z

    Arguments:
        z -- array of z variable
        vw -- wall velocities
        Dvw -- jump on the velocity profile at the wall
        sh -- determines if a shock is formed in the profile
        vsh -- value of the velocity at the - side of the shock
        xsh -- position of the shock

    Returns:
        fpzinf -- f'(z) function in the z -> infinity limit
        fpenv -- envelope of the f'(z) function in the z -> infinity limit
    '''

    fpzinf = -vw*np.sin(z*vw)*Dvw
    fpenv = vw*abs(Dvw)
    if sh:
        fpzinf += vsh*np.sin(z*xsh)*xsh
        fpenv += vsh*xsh

    return fpzinf*4*np.pi/z**2, fpenv*4*np.pi

########################### Constant profile ###########################
    
def fp_const(z, xib, xif, vconst=1):
    
    '''
    Function that computes the function f'(z) assuming a constant velocity
    profile with amplitude vconst and compact support between xib (back)
    and xif (front)
    
    Arguments:
        z -- array of z = k (t - t_n) values
        xib -- position of back shock
        xif -- position of front shock
        vconst -- amplitude of velocity profile (default is 1)

    Returns:
        fp_c -- analytical f'(z) for a constant profile
        fp_c0 -- value in the z -> 0 limit of f'(z)/z
    '''

    fp_c = xif*np.sin(z*xif) - xib*np.sin(z*xib)
    fp_c += 2/z*(np.cos(z*xif) - np.cos(z*xib))
    fp_c *= 4*np.pi*vconst/z**2
    fp_c0 = -np.pi*vconst/3*(xif**4 - xib**4)

    return fp_c, fp_c0

def vconst(xis, vs, n=1, xif=1, xib=0):
    
    '''
    Function that computes the characteristic constant velocity used for toy models
    from the n moment
    
    Note that n = 4 gives the amplitude of the constant profile that recovers the
    f'(z) value in the z -> 0 limit
    
    Arguments:
        xis -- array of xi values
        vs -- velocity profile
        n -- moment of the constant velocity (default is 1)
        xif -- position of the front shock (default is 1)
        xib -- position of the back shock (default is 0)
    '''

    vc = n*np.trapz(vs*xis**(n - 1), xis)/(xif**n - xib**n)

    return vc

def vchar(vw, vpl, vm, vsh=0, xish=0):
    
    '''
    Function that computes the characteristic constant velocity that leads to
    the same asymptotical limit in z to infinity of the f' function
    '''
    
    vc = (xish*vsh + vw*abs(vpl - vm))/(vw + xish)

    return vc

########################### Linear profiles ###########################

def fp_linear(z, xib, xif, v0=1, opt='inc'):
    
    '''
    Function that computes the function f'(z) assuming a linear velocity
    profile with maximum amplitude v0 and compact support between xib (back)
    and xif (front).
    The linear profile can either increase (opt='inc') or decrease linearly.
    '''
    
    Dxi = xif - xib

    # compute f'(z) for increasing linear profile
    if opt == 'inc':
        fp_l = xif*Dxi*np.sin(z*xif) - 1/z*(xib*np.cos(z*xib))
        fp_l += -1/z*(2*xib - 3*xif)*np.cos(z*xif)
        fp_l += -3/z**2*(np.sin(z*xif) - np.sin(z*xib))
        
    # compute f'(z) for decreasing linear profile
    if opt == 'dec':
        fp_l = -xib*Dxi*np.sin(z*xib) - 1/z*(xif*np.cos(z*xif))
        fp_l += 1/z*(3*xib - 2*xif)*np.cos(z*xib)
        fp_l += 3/z**2*(np.sin(z*xif) - np.sin(z*xib))
        
    fp_l *= 4*np.pi*v0/z**2/Dxi
    
    return fp_l

def vfp0_lin(xi, v, xif, xib, opt='inc'):
    
    '''
    Function that computes the max velocity of the linear profile such
    that f'(z) has the same value in the z -> 0 limit as the original profile
    '''
    
    # change xif by xib when the linear profile is decreasing
    if opt == 'inc':
        xif0 = xif
        xif = xib
        xib = xif0

    vfp0 = 20*(xif - xib)/(xif**5 + 4*xib**5 - 5*xif*xib**4)
    vfp0 *= np.trapz(xi**3*v, xi)
    
    return vfp0

######################### CORRELATION FUNCTIONS B_ij (r) #########################

#### this part is under development, to be checked and tested before made public

# compute bL and bN from f'(z)
def bn_bl_from_f(z, fp2, xi, multi=True, quiet=False):
    
    xi_ij, z_ij = np.meshgrid(xi, z, indexing='ij')
    zxi_ij = z_ij*xi_ij
    j1_z_oz = np.sin(zxi_ij)/zxi_ij**3 - np.cos(zxi_ij)/zxi_ij**2
    j0_z = np.sin(zxi_ij)/zxi_ij
    # avoid division by zero (use limits of functions for this value)
    j0_z[np.where(zxi_ij == 0)] = 1
    j1_z_oz[np.where(zxi_ij == 0)] = 1/3
    
    if multi:
        Nvws = np.shape(fp2)[0]
        BN = np.zeros((Nvws, len(xi)))
        BL = np.zeros((Nvws, len(xi)))
        
        for i in range(0, Nvws):
            xi_ij, fpz_ij = np.meshgrid(xi, fp2[i, :], indexing='ij')
            BN[i, :] = np.trapz(j1_z_oz*z_ij**2*fpz_ij**2, z, axis=1)/2/np.pi**2
            BL[i, :] = np.trapz((j0_z - 2*j1_z_oz)*z_ij**2*fpz_ij**2, z, axis=1)/2/np.pi**2
            if not quiet: print('vw ', i + 1, '/', Nvws, ' computed')

    else:
        xi_ij, fpz_ij = np.meshgrid(xi, fp2, indexing='ij')
        BN = np.trapz(j1_z_oz*z_ij**2*fpz_ij**2, z)/2/np.pi**2
        BL = np.trapz((j0_z - 2*j1_z_oz)*z_ij**2*fpz_ij**2, z, axis=1)/2/np.pi**2
        
    return BN, BL

# compute first and second derivatives of bL from f'(z)
def ders_bl_from_f(z, fp2, xi, multi=True, quiet=False):
    
    xi_ij, z_ij = np.meshgrid(xi, z, indexing='ij')
    zxi_ij = z_ij*xi_ij
    j0_z = np.sin(zxi_ij)/zxi_ij
    j1_z = np.sin(zxi_ij)/zxi_ij**2 - np.cos(zxi_ij)/zxi_ij
    
    integr1 = -j1_z - 2*j0_z/zxi_ij + 6*j1_z/zxi_ij**2
    integr2 = -j0_z + 4*j1_z/zxi_ij + 8*j0_z/zxi_ij**2 - 24*j1_z/zxi_ij**3
    # avoid division by zero (use limits of functions for this value)
    integr1[np.where(zxi_ij == 0)] = 0
    integr2[np.where(zxi_ij == 0)] = -1/5
    
    if multi:
        Nvws = np.shape(fp2)[0]
        der_BL = np.zeros((Nvws, len(xi)))
        der2_BL = np.zeros((Nvws, len(xi)))
        
        for i in range(0, Nvws):
            xi_ij, fpz_ij = np.meshgrid(xi, fp2[i, :], indexing='ij')
            der_BL[i, :] = np.trapz(integr1*z_ij**3*fpz_ij**2, z, axis=1)/2/np.pi**2
            der2_BL[i, :] = np.trapz(integr2*z_ij**4*fpz_ij**2, z, axis=1)/2/np.pi**2
            if not quiet: print('vw ', i + 1, '/', Nvws, ' computed')

    else:
        xi_ij, fpz_ij = np.meshgrid(xi, fp2, indexing='ij')
        der_BL[i, :] = np.trapz(integr1*z_ij**3*fpz_ij**2, z, axis=1)/2/np.pi**2
        der2_BL[i, :] = np.trapz(integr2*z_ij**4*fpz_ij**2, z, axis=1)/2/np.pi**2
        
    return der_BL, der2_BL

def correlations_from_f(xi, z, fp, first=True, second=True, quiet=False):
    
    '''
    Function to compute the lateral and normal component of the correlation
    function and their first and second derivative
    '''

    xi_ij, z_ij = np.meshgrid(xi, z, indexing='ij')
    zxi_ij = z_ij*xi_ij
    j0 = np.sin(zxi_ij)/zxi_ij
    j1 = np.sin(zxi_ij)/zxi_ij**2 - np.cos(zxi_ij)/zxi_ij
    integ0 = j0 - 2*j1/zxi_ij
    integN0 = j1/zxi_ij
    
    if second: first = True
    
    if first:
        integ_der0 = -j1 - 2*j0/zxi_ij + 6*j1/zxi_ij**2
    if second:
        integ_der20 = -j0 + 4*j1/zxi_ij + 8*j0/zxi_ij**2 - 24*j1/zxi_ij**3

    bl = np.zeros((len(vws), len(xi)))
    bn = np.zeros((len(vws), len(xi)))
    if first:
        der_bl = np.zeros((len(vws), len(xi)))
        der_bn = np.zeros((len(vws), len(xi)))
    if second:
        der2_bl = np.zeros((len(vws), len(xi)))
        der2_bn = np.zeros((len(vws), len(xi)))

    for i in range(0, len(vws)):

        xi_ij, fpz_ij = np.meshgrid(xi, fp[i, :], indexing='ij')
        integ = integ0*fpz_ij**2*z_ij**2
        integN = integN0*fpz_ij**2*z_ij**2
        if first: integ_der = integ_der0*fpz_ij**2*z_ij**3
        if second: integ_der2 = integ_der20*fpz_ij**2*z_ij**4

        bl[i, :] = np.trapz(integ, z, axis=1)/2/np.pi**2
        bn[i, :] = np.trapz(integN, z, axis=1)/2/np.pi**2
        if first:
            der_bl[i, :] = np.trapz(integ_der, z, axis=1)/2/np.pi**2
        if second:
            der2_bl[i, :] = np.trapz(integ_der2, z, axis=1)/2/np.pi**2
            
        if first:
            der_bn[i, :] = (bl[i, :] - bn[i, :])/xi
        if second:
            der2_bn[i, :] = (der_bl[i, :] - 2*der_bn[i, :])/xi
            
        if not quiet: print('vw ', i + 1, '/', Nvws, ' computed')
            
    if second:
        return bl, bn, der_bl, der_bn, der2_bl, der2_bn
    elif first:
        return bl, bn, der_bl, der_bn
    else:
        return bl, bn

def ders_bn_from_f(BL, BN, dBL, xi, multi=True, quiet=False):

    if multi:
        Nvws = np.shape(BL)[0]
        der_BN = np.zeros((Nvws, len(xi)))
        der2_BN = np.zeros((Nvws, len(xi)))
        
        for i in range(0, Nvws):
            der_BN[i, :] = (BL[i, :] - BN[i, :])/xi
            der2_BN[i, :] = (dBL[i, :] - 2*der_BN[i, :])/xi

    else:
        der_BN = (BL - BN)/xi
        der2_BN = (dBL - 2*der_BN)/xi

    return der_BN, der2_BN

def compute_bii_const(xir, xif, xib, retI=False):
    
    xip = xif + xib
    xid = xif - xib

    I1 = np.zeros(len(xir))
    I2 = np.zeros(len(xir))
    I3 = np.zeros(len(xir))

    for A in [-1, 1]:
        for B in [-1, 1]:

            x1 = xir - A*xif + B*xif
            I1 += xif**2*A*B*x1**2*np.sign(x1)
            I2 += xif*A*x1**3*np.sign(x1)
            I3 += x1**4*np.sign(x1)

            x2 = xir - A*xib + B*xib
            I1 += xib**2*A*B*x2**2*np.sign(x2)
            I2 += xib*A*x2**3*np.sign(x2)
            I3 += x2**4*np.sign(x2)

            x3 = xir - A*xif + B*xib
            I1 += -2*xif*xib*A*B*x3**2*np.sign(x3)
            I2 += -xif*A*x3**3*np.sign(x3)
            I3 += -2*x3**4*np.sign(x3)

            x4 = xir - A*xib + B*xif
            I2 += -xib*A*x4**3*np.sign(x4)

    I1 = -np.pi/16*I1
    I2 = np.pi/12*I2
    I3 = np.pi/48*I3
    
    bii = 8/xir*(I1 + I2 + I3)
    
    if retI:
        return bii, I1, I2, I3
    
    else:
        return bii