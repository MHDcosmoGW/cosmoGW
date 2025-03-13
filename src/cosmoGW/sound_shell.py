"""
sound_shell.py is a Python routine that contains functions used to compute the
kinetic spectrum produced from sound waves following the sound shell model.

Author: Alberto Roper Pol
Created: 01/02/2023

Other contributors: Antonino Midiri, Simona Procacci

Main reference is:

RPPC23 - A. Roper Pol, S. Procacci, C. Caprini,
"Characterization of the gravitational wave spectrum from sound waves within
the sound shell model," Phys. Rev. D 109, 063531 (2024), arXiv:2308.12943.
"""

import numpy as np
import matplotlib.pyplot as plt

import hydro_bubbles as hb
import plot_sets
import GW_analytical as an

### Reference values
cs2 = 1/3      # speed of sound
cols = hb.cols_ref

##### ROUTINES TO COMPUTE THE KINETIC SPECTRA USING THE SOUND-SHELL MODEL #####

def compute_kin_spec_dens(z, vws, fp, l, sp='sum', type_n='exp', cs2=cs2, min_qbeta=-4,
                          max_qbeta=5, Nqbeta=1000, min_TT=-1, max_TT=3, NTT=5000):
    
    '''
    Function that computes the kinetic power spectral density assuming exponential or simultaneous
    nucleation.
    
    Arguments:
        z -- array of values of z
        vws -- array of wall velocities
        fp -- function f'(z) computed from the sound-shell model using compute_fs
        l -- function lambda(z) computed from the sound-shell model using compute_fs
        sp -- type of function computed for the kinetic spectrum description
        type_n -- type of nucleation hystory (default is exponential 'exp',
                  another option is simultaneous 'sym')
        dens_spec -- option to return power spectral density (if True, default), or kinetic
                     spectrum (if False)
    '''
    
    if sp == 'sum': A2 = .25*(cs2*l**2 + fp**2)
    if sp == 'only_f': A2 = .5*fp**2               # assumes equal contributions from fp and cs l
    if sp == 'diff': A2 = .25*(fp**2 - cs2*l**2)
    if sp == 'cross': A2 = -.5*fp*np.sqrt(cs2)*l

    q_beta = np.logspace(min_qbeta, max_qbeta, Nqbeta)
    TT = np.logspace(min_TT, max_TT, NTT)
    q_ij, TT_ij = np.meshgrid(q_beta, TT, indexing='ij')
    Pv = np.zeros((len(vws), len(q_beta)))

    funcT = np.zeros((len(vws), len(q_beta), len(TT)))
    for i in range(0, len(vws)):
        if type_n == 'exp':
            funcT[i, :, :] = np.exp(-TT_ij)*TT_ij**6*np.interp(TT_ij*q_ij, z, A2[i, :])
        if type_n == 'sim':
            funcT[i, :, :] = .5*np.exp(-TT_ij**3/6)*TT_ij**8*np.interp(TT_ij*q_ij, z, A2[i, :])
        Pv[i, :] = np.trapz(funcT[i, :, :], TT, axis=1)

    return q_beta, Pv

# def compute_kin_spec(vws, q_beta, Pv, corr=True, cs2=cs2):
    
#     '''
#     Function that computes the kinetic spectrum as a function of k Rast from
#     the power spectral density as a function of q/beta
#     '''

#     EK = np.zeros((len(vws), len(q_beta)))
#     cs = np.sqrt(cs2)
#     if corr: Rstar_beta = (8*np.pi)**(1/3)*np.maximum(vws, cs)
#     else: Rstar_beta = (8*np.pi)**(1/3)*vws
#     kks = np.zeros((len(vws), len(q_beta)))

#     for i in range(0, len(vws)):
#         kks[i, :] = q_beta*Rstar_beta[i]
#         pref = kks[i, :]**2/Rstar_beta[i]**6/(2*np.pi**2)
#         EK[i, :] = pref*Pv[i, :]

#     return kks, EK

def compute_kin_spec(vws, q_beta, Pv, corr=True, cs2=cs2):

    EK = []
    cs = np.sqrt(cs2)
    if corr: Rstar_beta = (8*np.pi)**(1/3)*np.maximum(vws, cs)
    else: Rstar_beta = (8*np.pi)**(1/3)*vws
    kks = np.zeros((len(vws), len(q_beta)))

    for i in range(0, len(vws)):
        kks[i, :] = q_beta*Rstar_beta[i]
        pref = kks[i, :]**2/Rstar_beta[i]**6/(2*np.pi**2)
        EK_aux = pref*Pv[i, :]
        EK_aux2 = kin_sp(kks[i, :], EK_aux, cs2=cs2)
        EK = np.append(EK, EK_aux2)

    return kks, EK

########################### NOT PUBLIC ON GITHUB ###########################

def compute_lamb_vip(z, xi, v, lamb, cs2=cs2):

    xi_ij, z_ij = np.meshgrid(xi, z, indexing='ij')
    lamb_ij, z_ij = np.meshgrid(lamb, z, indexing='ij')
    v_ij, z_ij = np.meshgrid(v, z, indexing='ij')

    l_z = 4*np.pi/z*np.trapz(lamb_ij*xi_ij*np.sin(xi_ij*z_ij), xi, axis=0)
    lp_z = -l_z/z + 4*np.pi/z*np.trapz(lamb_ij*xi_ij**2*np.cos(xi_ij*z_ij), xi, axis=0)
    
    f_z = 4*np.pi/z*np.trapz(v_ij*np.sin(z_ij*xi_ij), xi, axis=0)
    fp_z = -f_z/z + 4*np.pi/z*np.trapz(xi_ij*v_ij*np.cos(z_ij*xi_ij), xi, axis=0)
    
    gp_z = fp_z - 4*np.pi*np.trapz(xi_ij**2*v_ij*np.sin(z_ij*xi_ij), xi, axis=0)
    
    return f_z, fp_z, gp_z, l_z, lp_z

class funcs():
    
    def __init__(self, z, f, l, fp, gp, lp, vws, alp):
        self.z = z
        self.f = f
        self.l = l
        self.fp = fp
        self.gp = gp
        self.lp = lp
        self.vws = vws
        self.alpha = alp

def compute_fs(vws, alpha, bag_eos=True, cs2=cs2, Nxi=10000, Nz=10000, z_min=-3, z_max=5, vvs_bey_beos=[],
               wws_bey_beos=[], save=True, alp_str='', read=False):
    
    '''
    Function that computes the object ff, which contains relevant functions used
    in the sound shell model: f(z), g(z), lambda(z) and derivatives from the
    integrals over the 1d profiles of the velocity and enthalpy of expanding 
    bubbles.
    
    Arguments:
        vws -- array of wall velocities
        alpha -- value of alpha at the nuccleation temperature
        bag_eos -- option to choose the bag eos (at the moment, this is the only option
                   available for computation, if bag_eos is chosen to be False, then one
                   needs to provide 1d profiles in vvs_bey_beos, wws_bey_beos)
        cs2 -- square of the speed of sound (default is 1/3)
        Nxi -- number of points in xi discretization
        Nz -- number of points in z discretization
        z_min -- minimum value of z in log (default is 1e-3)
        z_max -- maximum value of z in log (default 1e5)
        save -- option to save the results in a pickle variable under 'results/alpha_%alpha/ff.pckl'
        read -- option to directly read the results from previously saved pickle variable
    '''
    
#    if bag_eos:
    
    #xis, vvs, lams, alphas_pl, conv = hb.compute_profiles_vws(alpha, cs2=cs2, vws=vws,
    #                                               Nxi=Nxi, plot=False, alphan=True, plot_v='both',
    #                                               tol=1e-4, quiet=True, lam=True, legs=True, st_lg=1)
    
    if alp_str == '': alp_str = str(alpha)
    
    if not read:

        xis, vvs, wws, alphas_pl, conv = hb.compute_profiles_vws(alpha, cs2=cs2, vws=vws,
                                                       Nxi=Nxi, plot=False, alphan=True, plot_v='both',
                                                       tol=1e-4, quiet=True, lam=False, legs=True, st_lg=1)

    # if one uses cs^2 different than 1/3 (beyond bag eos), then profiles can be given as an input
    # to this function
    if not bag_eos:
        vvs = vvs_bey_beos
        wws = wws_bey_beos
    
    if not read:
        # compute functions f', g', l, f
        z = np.logspace(z_min, z_max, Nz)

        fzs = np.zeros((len(vws), len(z)))
        fpzs = np.zeros((len(vws), len(z)))
        gpzs = np.zeros((len(vws), len(z)))
        lzs = np.zeros((len(vws), len(z)))
        lpzs = np.zeros((len(vws), len(z)))
        lams = np.zeros((len(vws), len(xis)))
        kappas = np.zeros(len(vws))

        for ind in range(0, len(vws)):
            lams[ind, :] = hb.w_to_lam(xis, wws[ind, :], vws[ind], alpha)
            fz, fpz, gpz, lz, lpz = compute_lamb_vip(z, xis, vvs[ind, :], lams[ind, :])
            fzs[ind, :] = fz
            fpzs[ind, :] = fpz
            gpzs[ind, :] = gpz
            lzs[ind, :] = lz
            lpzs[ind, :] = lpz
            kappas[ind] = 4/vws[ind]**3/alpha*np.trapz(xis**2*wws[ind, :]/(1 - vvs[ind, :]**2)*vvs[ind, :]**2, xis)

            print('Computed ', ind, '/', len(vws) - 1)

        ff = funcs(z, fzs, lzs, fpzs, gpzs, lpzs, vws, alpha)
        ff.kappa = kappas

    else:

        import pickle
        fll = 'results/alpha_%s/ff.pckl'%alp_str
        with open(fll, 'rb') as handle:
            ff = pickle.load(handle)
        print('Read pickle variable containing ff in ', fll)


    if save and not read:
        import pickle
        fll = 'results/alpha_%s/ff.pckl'%alp_str
        with open(fll, 'wb') as handle:
            pickle.dump(ff, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print('Saved pickle variable containing ff in ', fll)

    return ff

def compute_fp_const(z, xib, xif):
    
    fp = 4*np.pi/z**2*(xif*np.sin(z*xif) - xib*np.sin(z*xib) + 2/z*(np.cos(z*xif) - np.cos(z*xib)))
    
    return fp

### Function added by Antonino for computing kinetic spectra before collisions

def compute_kin_spec_densBC_NoDeltaApprox(z, vws, fp, l, sp='sum', type_n='sim', limit='full', cs2=cs2, min_qbeta=-4,
                          max_qbeta=5, Nqbeta=1000, tcrit=-10, tnucl=0, tcurr=10, NTT=5000, beta_2=1, beta_eff=1, normexp=1, normsim=1):
    
    '''
    Function that computes the kinetic power spectral density assuming exponential or simultaneous
    nucleation BEFORE COLLISIONS (AVERAGING OVER NUCLEATION TIMES) -- Here the simultaneous case is treated without taking the
    large time limit in which we can use the delta approximation
    
    Arguments:
        z -- array of values of z
        vws -- array of wall velocities
        fp -- function f'(z) computed from the sound-shell model using compute_fs
        l -- actually not needed here
        sp -- type of function computed for the kinetic spectrum description
        type_n -- type of nucleation hystory (default is exponential 'exp',
                  another option is simultaneous 'sym')
        dens_spec -- option to return power spectral density (if True, default), or kinetic
                     spectrum (if False)
        tcrit -- critical time
        tnucl -- nucleation time (or t_* in the draft considering that p(t) is expanded around this time)
        tcurr -- current time
        beta_2 -- parameter in the draft
        beta_eff -- parameter in the draft
        normexp, normsim -- ''-ln(h(t_*))'' (can be changed according to the desired normalization)
    '''
    
    cs = np.sqrt(cs2)
    
    if sp == 'sum': A2 = .25*(cs2*l**2 + fp**2)
    if sp == 'only_f': A2 = .25*fp**2
    if sp == 'diff': A2 = .25*(fp**2 - cs2*l**2)
    if sp == 'cross': A2 = .5*fp*cs*l

    q_beta = np.logspace(min_qbeta, max_qbeta, Nqbeta)
    TT = np.linspace(tcrit, tcurr, NTT)
    TTn = np.repeat(tnucl, NTT)
    TTc = np.repeat(tcurr, NTT)
    intint = np.linspace(beta_2*(tcrit-tnucl), beta_2*(tcurr-tnucl), 10000)
    q_ij, TT_ij = np.meshgrid(q_beta, TT, indexing='ij')
    q_ij, TTn_ij = np.meshgrid(q_beta, TTn, indexing='ij') 
    q_ij, TTc_ij = np.meshgrid(q_beta, TTc, indexing='ij')
    Pv = np.zeros((len(vws), len(q_beta)))

    funcT = np.zeros((len(vws), len(q_beta), len(TT)))
    for i in range(0, len(vws)):
        if type_n == 'exp':
            funcT[i, :, :] = np.exp(-normexp*np.exp(TT_ij-TTn_ij))*normexp*np.exp(TT_ij-TTn_ij)*(TTc_ij-TT_ij)**6*np.interp((TTc_ij-TT_ij)*q_ij, z, A2[i, :])
            Pv[i, :] = np.trapz(funcT[i, :, :], TT, axis=1)
        if type_n == 'sim':
            #Integral does the computation in the asymptotic limit in which we have exp(-beta_eff^3 * (t-t_n)^3 / 6)
            def Integral(i, beta_2, beta_eff):
                return np.trapz(normsim*np.exp(-beta_2**2*(TT_ij-TTn_ij)**2/2)*np.exp(-beta_eff**3*(TT_ij-TTn_ij)**3/6)*(TTc_ij-TT_ij)**6*np.interp((TTc_ij-TT_ij)*q_ij, z, A2[i, :])/np.sqrt(2*np.pi), TT, axis=1)
            #IntegralFull does the full computation of h(t) without the large time approximation
            def IntegralFull(i, beta_2, beta_eff):
                return np.trapz(normsim*np.exp(-beta_2**2*(TT_ij-TTn_ij)**2/2)*np.exp(-beta_eff**3/(6*np.sqrt(2*np.pi)*beta_2**3)*np.trapz(np.exp(-0.5*intint**2)*(beta_2*(TT_ij-TTn_ij)-intint)**3, intint))*(TTc_ij-TT_ij)**6*np.interp((TTc_ij-TT_ij)*q_ij, z, A2[i, :])/np.sqrt(2*np.pi), TT, axis=1)

            if limit == 'asymptotic':
                Pv[i, :] = Integral(i, beta_2, beta_eff)
            if limit == 'full':
                Pv[i, :] = IntegralFull(i, beta_2, beta_eff)
    
    return q_beta, Pv

### Function added by Antonino for computing kinetic spectra before collisions

def compute_kin_spec_densBC_GaussMollifier(z, vws, fp, l, sp='sum', type_n='sim', cs2=cs2, min_qbeta=-4,
                          max_qbeta=5, Nqbeta=1000, tcrit=-10, tnucl=0, tcurr=10, NTT=5000, sigma_ext=1):
    
    '''
    Function that computes the kinetic power spectral density assuming exponential or simultaneous
    nucleation BEFORE COLLISIONS (AVERAGING OVER NUCLEATION TIMES) -- Here the simultaneous case is treated by substituting
    the final delta function by a normalized gaussian of desired sigma (sigma_ext parameter)
    -- Not useful for the purpose of the paper (just a check)
    
    Arguments:
        z -- array of values of z
        vws -- array of wall velocities
        fp -- function f'(z) computed from the sound-shell model using compute_fs
        l -- function lambda(z) computed from the sound-shell model using compute_fs
        sp -- type of function computed for the kinetic spectrum description
        type_n -- type of nucleation hystory (default is exponential 'exp',
                  another option is simultaneous 'sym')
        dens_spec -- option to return power spectral density (if True, default), or kinetic
                     spectrum (if False)
        sigma_ext -- gives the sigma of the Gaussian by which we approximate the delta function in the simultaneous case (approaching the 
                     delta in the limit sigma_ext -> 0)
    '''
    
    cs = np.sqrt(cs2)
    
    if sp == 'sum': A2 = .25*(cs2*l**2 + fp**2)
    if sp == 'only_f': A2 = .25*fp**2
    if sp == 'diff': A2 = .25*(fp**2 - cs2*l**2)
    if sp == 'cross': A2 = .5*fp*cs*l

    q_beta = np.logspace(min_qbeta, max_qbeta, Nqbeta)
    TT = np.linspace(tcrit, tcurr, NTT)
    TTn = np.repeat(tnucl, NTT)
    TTc = np.repeat(tcurr, NTT)
    q_ij, TT_ij = np.meshgrid(q_beta, TT, indexing='ij')
    q_ij, TTn_ij = np.meshgrid(q_beta, TTn, indexing='ij') 
    q_ij, TTc_ij = np.meshgrid(q_beta, TTc, indexing='ij')
    Pv = np.zeros((len(vws), len(q_beta)))

    funcT = np.zeros((len(vws), len(q_beta), len(TT)))
    for i in range(0, len(vws)):
        if type_n == 'exp':
            funcT[i, :, :] = np.exp(-np.exp(TT_ij-TTn_ij))*np.exp(TT_ij-TTn_ij)*(TTc_ij-TT_ij)**6*np.interp((TTc_ij-TT_ij)*q_ij, z, A2[i, :])
            Pv[i, :] = np.trapz(funcT[i, :, :], TT, axis=1)
        if type_n == 'sim':
            def Integral(i, sigma):
                return np.trapz(1/(np.sqrt(2*np.pi)*sigma)*np.exp(-(TT_ij-TTn_ij)**2/(2*sigma**2))*(TTc_ij-TT_ij)**6*np.interp((TTc_ij-TT_ij)*q_ij, z, A2[i, :]), TT, axis=1)
            Pv[i, :] = Integral(i, sigma_ext)
    
    return q_beta, Pv

class kin_sp():
    
    def __init__(self, k, EK, cs2=cs2):
        self.EK = EK
        self.k = k
        ind_max = np.argmax(EK)
        self.max_EK = EK[ind_max]
        self.kmax = k[ind_max]
        self.zeta = EK/self.max_EK
        self.KK = np.trapz(self.zeta, k)
        self.CC = 16/15*np.trapz(self.zeta**2/k**2, k)
        self.CC2 = 16/15*np.trapz(self.zeta**2/k**4, k)
        self.DD = 16/15*np.trapz(self.zeta*k**2, k)/self.CC
        self.OmK = self.KK*self.max_EK
        self.amplGW = 3*(1 + cs2)*self.OmK**2/self.KK**2*self.CC

################## ROUTINES TO COMPUTE THE STRESSES SPECTRA ##################

def compute_zeta_Pi(kk, zetaK, kps, zetaPi, vws=[], inds_vws=[], Np=3000, Nk=120, plot=False,
                    largek=3, smallk=-3, extend=False, alp_str='01', typ='exp',
                    dir0='results', read=False, save=True, opt_stt=0):
    
    '''
    Function that computes the normalized anisotropic stresses spectrum for the kinetic
    spectra computed within the sound-shell model for fully compressional modes.
    
    It uses the function GW_analytical.ET_correlator_compr()
    
    Reference: A. Roper Pol, S. Procacci, C. Caprini, "Characterization of the
    gravitational wave spectrum from sound waves within the sound shell model,"
    https://arxiv.org/abs/2308.12943, eq. 47
    '''
    
    import pickle
    
    if inds_vws == []: inds_vws = range(0, len(vws))
    
    for i in inds_vws:
        
        dirr = dir0 + '/alpha_%s/'%alp_str

        if opt_stt < 1: stt = '%.0f'%(10*np.around(vws[i], decimals=1))
        else: stt = '%.0f'%(100*np.around(vws[i], decimals=2))
        if stt != '10' and opt_stt < 1:
            stt = '0' + stt
        strr = dirr + 'zetaPi_%s_%s'%(typ, stt)
        strk = dirr + 'kps_%s'%stt
        if stt == '10':
            strr = dirr + 'zetaPi_%s_1'%typ
            strk = dirr + 'kps_1'
        
        if read:

            with open(strk, 'rb') as handle:
                kps[i, :] = pickle.load(handle)
                
            print('read k array from ', strk)

            with open(strr, 'rb') as handle:
                zetaPi[i, :] = pickle.load(handle)
                
            print('read zeta_Pi array from ', strr)
            
        else:

            kps[i, :], zetaPi[i, :], CC = \
                    an.ET_correlator_compr(kk[i, :], zetaK[i, :], Np=Np, Nk=Nk, plot=plot,
                                           extend=extend, largek=largek, smallk=smallk)
            
            print('computed stress ', i + 1, '/', len(inds_vws))
            
            if save:

                with open(strk, 'wb') as handle:
                    pickle.dump(kps[i, :], handle, protocol=pickle.HIGHEST_PROTOCOL)
                    
                print('saving k array in ', strk)
                    
                with open(strr, 'wb') as handle:
                    pickle.dump(zetaPi[i, :], handle, protocol=pickle.HIGHEST_PROTOCOL)
                    
                print('saving zeta_Pi array in ', strr)

    return kps, zetaPi

def plot_zeta_Pi(kps, zeta_Pi, vws=[], inds_vws=[], compk=0, typ='exp',
                 x0=1e-2, x1=1e3, y0=1e-9, y1=4, locc='upper left', save=True,
                 dir0='plots', opt_stt=0):
    
    if inds_vws == []: inds_vws = range(0, len(vws))
    
    kGW = np.zeros(len(vws))
    k3Pi_max = np.zeros(len(vws))

    for ind in inds_vws:
        
        ind_c = ind%11
        
        lbl = '%.1f'%vws[ind]
        if opt_stt >= 1: lbl = '%.2f'%vws[ind]

        plt.plot(kps[ind, :], kps[ind, :]**compk*zeta_Pi[ind, :],
                 color=cols[ind_c], label=r'$\xi_w = %s$'%lbl)
    
        good = np.where(kps[ind, :] < 1e3)[0]
        ind_max = np.argmax(kps[ind, good]**3*zeta_Pi[ind, good])
        kGW[ind] = kps[ind, good][ind_max]
        k3Pi_max[ind] = kGW[ind]**compk*zeta_Pi[ind, good][ind_max]
        plt.plot(kGW[ind], k3Pi_max[ind], color=cols[ind_c], marker='*',
                 markersize=10)

    plt.loglog()
    plt.legend(frameon=False, fontsize=13, loc=locc)

    plt.xlabel(r'$K \equiv kR_*$')
    pp = ''
    if compk != 0: pp = 'K^%i'%compk
    plt.ylabel(r'$%s \zeta_{\Pi} (K)$'%pp)
    plot_sets.axes_lines()

    plt.xlim(x0, x1)
    plt.ylim(y0, y1)
    
    if save:
        fll = dir0 + '/'
        if compk != 0: fll += 'K%i_'%compk
        fll += 'zeta_Pi_%s.pdf'%typ
        plt.savefig(fll, bbox_inches='tight')
        print('saved plot in ', fll)
    
    return kGW, k3Pi_max

################## ROUTINES TO COMPUTE THE GW SPECTRA ##################

def compute_zeta_GW(kk, zetaK, kps, zetaGW, vws=[], inds_vws=[], kHs=[], inds_kHs=[],
                    dtfs=[], inds_dtfs=[], Np=3000, Nk=120, plot=False, largek=3, smallk=-3,
                    expansion=True, cs2=cs2, terms='all', inds_m=[], inds_n=[], corr_Delta_0=True,
                    tini=1., extend=False, alp_str='01', typ='exp', dir0='results', read=False, save=True,
                    cs2_str='', quiet=False, opt_stt=0):

    '''
    Function that computes the normalized GW spectrum for the kinetic
    spectra computed within the sound-shell model for fully compressional modes
    assuming stationary UETC.
    
    It uses the function GW_analytical.effective_ET_correlator_stat()
    
    Reference: A. Roper Pol, S. Procacci, C. Caprini, "Characterization of the
    gravitational wave spectrum from sound waves within the sound shell model,"
    https://arxiv.org/abs/2308.12943, eq. 90
    '''
    
    import pickle
    
    if inds_vws == []: inds_vws = range(0, len(vws))
    if inds_dtfs == []: inds_dtfs = range(0, len(dtfs))
    if inds_kHs == []: inds_kHs = range(0, len(kHs))
    
    dirr = dir0 + '/alpha_%s/'%alp_str

    for i in inds_vws:
        if not quiet: print('vw: ', vws[i], ', ', i + 1, '/', len(vws))
        for j in inds_dtfs:
            if not quiet: print('dtfs: ', dtfs[j], ', ', j + 1, '/', len(dtfs))
            for m in inds_kHs:
                if not quiet: print('kH: ', kHs[m], ', ', m + 1, '/', len(kHs))
                
                if opt_stt < 1: stt = '%.0f'%(10*np.around(vws[i], decimals=1))
                else: stt = '%.0f'%(100*np.around(vws[i], decimals=2))
                if stt != '10' and opt_stt < 1:
                    stt = '0' + stt
                strr = dirr + 'zetaGW_%s_%s'%(typ, stt)
                if stt == '10': strr = dirr + 'zetaGW_%s_1'%typ
                
                strr += '_kH%i'%kHs[m]
                
                exp = np.floor(np.log10(dtfs[j]))
                if opt_stt < 1: strr += '_dtf1e%i'%np.log10(dtfs[j])
                else: strr += '_dtf%.1fe%i'%(dtfs[j]/10**exp, exp)
                
                if cs2_str != '': strr += '_cs2_' + cs2_str
                
                if read:

                    with open(strr, 'rb') as handle:
                        zetaGW[i, j, m, :] = pickle.load(handle)
                        
                    if not quiet: print('read GW spectrum from ', strr)
            
                else:

                    kps[i, :], zetaGW[i, j, m, :], CC = \
                            an.effective_ET_correlator_stat(kk[i, :], zetaK[i, :], dtfs[j] + tini, Np=Np, Nk=Nk,
                                                            plot=plot, expansion=expansion, terms=terms,
                                                            inds_m=inds_m, inds_n=inds_n, kstar=kHs[m], extend=extend,
                                                            largek=largek, smallk=smallk, tini=tini, cs2=cs2,
                                                            corr_Delta_0=corr_Delta_0)

                    if not quiet: print('computed GW spectrum for kHs = ', kHs[m], ', dtfs = ', dtfs[j],
                          ', vws = ', vws[i])
            
                    if save:

                        with open(strr, 'wb') as handle:
                            pickle.dump(zetaGW[i, j, m, :], handle, protocol=pickle.HIGHEST_PROTOCOL)

                        if not quiet: print('saving zeta_GW array in ', strr)

    return kps, zetaGW