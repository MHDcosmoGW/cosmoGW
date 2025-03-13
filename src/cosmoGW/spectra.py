"""
spectra.py is a Python routine that contains description for specific spectral templates,
postprocessing routines for numerical spectra, and other mathematical routines.

Author: Alberto Roper Pol
Created: 01/01/2021
Updated: 03/09/2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

##### FUNCTIONS FOR SMOOTHED DOUBLE BROKEN POWER LAW SPECTRA #####

def smoothed_dbPL(k, a=4, b=5/3, alp=2, kpeak=1, Omega=False):
    
    """
    Function that returns the value of the smoothed double broken power law (dbPL) model
    for the magnetic or kinetic spectrum, which has the form
    
        zeta(K) = (b + abs(a))^(1/alp) K^a/[ b + c K^(alp(a + b)) ]^(1/alp),
        
    where K = k/kpeak, c = 1 if a = 0 or c = abs(a) otherwise.
    
    The function is only correctly defined when b > 0 and a + b >= 0.

    Arguments:
        x -- values of x
        a -- slope of the spectrum at low wave numbers, k^a
        b -- slope of the spectrum at high wave numbers, k^(-b)
        alp -- smoothness of the transition from one power law to the other
        kpeak -- spectral peak, i.e., position of the break from k^a to k^(-b)
        Omega -- option to use the integrated energy desnity as the input A

    Returns:
        y -- value of the spectrum
    """
    
    if b < max(0, -a):
        print('b has to be larger than 0 and -a')
        return 0*k**0
    
    c = abs(a)
    if a == 0: c = 1
    m = (b + abs(a))**(1/alp)
    
    return m*(k/kpeak)**a/(b + c*(k/kpeak)**((a+b)*alp))**(1/alp)

## The complete beta function is used to compute moments of the 
## smoothed dbPL spectra:
##
## \int zeta(K) K^n dK
##
## The integral converges only when a + n > -1 and b - n > 1
##

def complete_beta(a, b):
    
    import math as m
    return m.gamma(a)*m.gamma(b)/m.gamma(a + b)

def calIab_n_alpha(alp=2, a=4, b=5/3, n=0):
    
    '''
    Function that computes the prefactor of the integral
    
    Iab_n (alpha) = \int dK K^n zeta(K),
    
    where zeta(K) is the function defined in smoothed_dbPL in terms
    of incomplete beta functions, used in Iab_n_alpha and Iab_n_alpha_inc.
    '''

    c = abs(a)
    if a == 0: c = 1

    alp2 = 1/alp/(a + b)
    Iab = ((abs(a) + b)/b)**(1/alp)*alp2
    Iab *= (b/c)**((a + n + 1)*alp2)

    return Iab

def Iabn(alp=2, a=4, b=5/3, n=0):
    
    '''
    Function that computes the moment n of the smoothed dbPL spectra:
    
    (b + abs(a))^(1/alp) \int K^(a + n)/(b + c K^(alp(a + b))) dK,
    
    where c = 1 if a = 0 or c = abs(a) otherwise.
    '''
    
    alp2 = 1/alp/(a + b)
    a_beta = (a + n + 1)*alp2
    b_beta = (b - n - 1)*alp2
    
    if b_beta > 0 and a_beta > 0:
            
        calI = calIab_n_alpha(alp=alp, a=a, b=b, n=n)
        comp_beta = complete_beta(a_beta, b_beta)
        return comp_beta*calI
    
    else:
        if b_beta <= 0:
            print('b + n has to be larger than 1 for the integral',
                  'to converge')
        if a_beta <= 0:
            print('a + n has to be larger than -1 for the integral',
                  'to converge')
        return 0

def calA(alp=2, a=4, b=5/3):
    
    '''
    Function that computes the parameter {\cal A} = Iab,0 that relates the
    peak and the integrated values of the smoothed_dbPL spectrum
    '''
    
    return Iabn(alp=alp, a=a, b=b, n=0)

def calB(alp=2, a=4, b=5/3):
    
    '''
    Function that computes the parameter {\cal B} = Iab,-1/Iab,0 that relates the
    peak and the integral scale
    '''
    
    BB = Iabn(alp=alp, a=a, b=b, n=-1)/Iabn(alp=alp, a=a, b=b, n=0)
    
    return BB

def Iabn_infty(a=4, b=5/3, n=0):
    
    '''
    Function that computes the moment n of the piece-wise broken power law
    zeta(K) = K^a if K <= 1 and K^{-b} if K > 1
    '''
    
    return (a + b)/(a + n + 1)/(b - n - 1)

def calB_infty(a=4, b=5/3):
    
    '''
    Function that computes the the parameter {\cal B} = Iab,-1/Iab,0 that relates the
    peak and the integral scale for the piece-wise broken power law
    
    zeta(K) = K^a if K <= 1 and K^{-b} if K > 1
    
    '''
    
    return (a + 1)*(b - 1)/a/b

def calB_general(alp=2, a=4, b=5/3, n=1):
    
    '''
    Function that computes the parameter {\cal B}_{ab,n} = (Iab,-n/Iab,0)^(1/n) that relates the
    peak and the integral scale of the n-th order
    '''
    
    BB = (Iabn(alp=alp, a=a, b=b, n=-n)/Iabn(alp=alp, a=a, b=b, n=0))**(1/n)
    
    return BB

def peak_comp_sp_dbpl(a=4, b=5/3, alp=2, n=1):
    
    '''
    Function that computes the position and the maximum of the zeta(P)
    smoothed double broken power law multiplied by k^n
    '''
    
    c = abs(a)
    if a == 0: c = 1
    
    kpos = (b/c*(a + n)/(b - n))**(1/alp/(a + b))
    Om_peak = ((abs(a) + b)*(b - n)/b/(a + b))**(1/alp)
    Om_peak *= (b*(a + n)/c/(b - n))**((a + n)/alp/(a + b))
    
    return kpos, Om_peak

def calC(alp=2, a=4, b=5/3):
    
    '''
    Function that computes the parameter {\cal C} = 28/15 * I2a2b,-2 that allows to
    compute the stress spectrum by taking the convolution of the smoothed_dbPL
    spectra over \kk and \kk - \pp.
    
    It gives the spectrum of the stress of Gaussian vortical non-helical fields
    as
    
    P_\Pi (0) = 2 \pi^2 EM*^2 {\cal C} / k*
    
    '''
    
    return 28/15*Iabn(alp=alp/2, a=a*2, b=b*2, n=-2)

def calCinf(alp=2, a=4, b=5/3):
    
    '''
    Function that computes the parameter {\cal C}^\infty that allows to
    compute the asymptotic branch of the stress spectrum of the smoothed_dbPL
    for Gaussian vortical non-helical fields:
    
    p_\Pi (\infty) = {\cal C}^\infty K^{-b - 2}
    
    '''
    
    C = calC(alp=alp, a=a, b=b)
    A = calA(alp=alp, a=a, b=b)
    
    return 8/3*((b + abs(a))/abs(a))**(1/alp)*A/C

def calC_infty(a=4, b=5/3):
    
    '''
    Function that computes the parameter {\cal C} = 28/15 * I2a2b,-2 that allows to
    compute the stress spectrum by taking the convolution of the smoothed_dbPL
    spectra over \kk and \kk - \pp.
    
    It gives the spectrum of the stress of Gaussian vortical non-helical fields
    as
    
    P_\Pi (0) = 2 \pi^2 EM*^2 {\cal C} / k*
    
    '''
    
    return 56/15*(a + b)/(2*a - 1)/(2*b + 1)

def calCinf_infty(a=4, b=5/3):
    
    '''
    Function that computes the parameter {\cal C} = 28/15 * I2a2b,-2 that allows to
    compute the stress spectrum by taking the convolution of the smoothed_dbPL
    spectra over \kk and \kk - \pp.
    
    It gives the spectrum of the stress of Gaussian vortical non-helical fields
    as
    
    P_\Pi (0) = 2 \pi^2 EM*^2 {\cal C} / k*
    
    '''
    
    return 5/7*(2*a - 1)*(2*b + 1)/(a + 1)/(b - 1)

######################### GENERAL FUNCTIONS ################################

def compute_kpeak(k, E, tol=.01, quiet=False):

    """
    Function that computes the maximum of the spectrum E and its spectral
    peak.

    Arguments:
        k -- array of wave numbers
        E -- array of the spectral values
        tol -- factor to avoid faulty maxima due to nearly flat spectrum
               (default 1%)
        quiet -- option to print out the result if quiet is False
                 (default False)

    Return:
        kpeak -- position of the spectral peak
        Emax -- maximum value of the spectrum
    """

    max1 = np.argmax(E)
    indmax = max1

    if E[max1] == 0:
        Emax = 0
        kpeak = 0
    else:
        max2 = np.argmax(k*E)
        # if the maximum of the spectrum is within tol of the maximum of k*E,
        # then we take as the maximum value where k*E is maximum, to take into
        # account flat and nearly flat spectra
        if abs(E[max1] - E[max2])/E[max1] < tol: indmax = max2
        Emax = E[indmax]
        kpeak = k[indmax]

    if not quiet:
        print('The maximum value of the spectrum is ', Emax,
              ' and the spectral peak is ', kpeak)

    return kpeak, Emax

def max_E_kf(k, E, exp=0):

    """
    Function that computes the maximum of a spectrum compensated by the
    wave number, i.e., max(E*k^exp)

    Arguments:
        k -- array of wave numbers
        E -- array of spectral values
        exp -- exponent of k (default 0)
    """

    indmax = np.argmax(k**exp*E)
    max_k = k[indmax]
    max_E = E[indmax]

    return max_k, max_E

def characteristic_k(k, E, exp=1):

    """
    Function that computes the characteristic wave number.

    Arguments:
        k -- array of wave numbers
        E -- array of spectral values
        exp -- exponent used to define the characteristic wave number
               k_ch ~ (\int k^exp E dk/\int E dk)^(1/exp)
               (default 1)

    Returns:
        kch -- characteristic wave number defined with the power 'exp'
    """

    k = k[np.where(k != 0)]
    E = abs(E[np.where(k != 0)])
    spec_mean = np.trapz(E, k)
    int = np.trapz(E*k**exp, k)
    # avoid zero division
    if exp >= 0 and spec_mean == 0: spec_mean = 1e-30
    if exp < 0 and int == 0: int = 1e-30
    kch = (int/spec_mean)**(1/exp)

    return kch

def min_max_stat(t, k, E, abs_b=True, indt=0, plot=False, hel=False):

    """
    Function that computes the minimum, the maximum, and the averaged
    functions over time of a spectral function.

    Arguments:
        t -- time array
        k -- wave number array
        E -- spectrum 2d array (first index t, second index k)
        indt -- index of time array to perform the average
                from t[indt] to t[-1]
        plot -- option to overplot all spectral functions
                from t[indt] to t[-1]
        hel -- option for helical spectral functions where positive and
               negative values can appear (default False)
               It then returns min_E_pos, min_E_neg, max_E_pos, max_E_neg
               referring to the maximum/minimum absolute values of the positive
               and negative values of the helical funtion.

    Returns:
        min_E -- maximum values of the spectral function over time
        max_E -- maximum values of the spectral function over time
        stat_E -- averaged values of the spectral function over time
                   from t[indt] to t[-1]
    """

    import matplotlib.pyplot as plt

    if hel:
        min_E_neg = np.zeros(len(k)) + 1e30
        max_E_neg = np.zeros(len(k))
        min_E_pos = np.zeros(len(k)) + 1e30
        max_E_pos = np.zeros(len(k))
    else:
        min_E = np.zeros(len(k)) + 1e30
        max_E = np.zeros(len(k))
    for i in range(indt, len(t)):
        if hel:
            if plot: plt.plot(k, abs(E[i,:]))
            # split between positive and negative values
            x_pos, x_neg, f_pos, f_neg, color = red_blue_func(k, E[i, :])
            for j in range(0, len(k)):
                if k[j] in x_pos:
                    indx = np.where(x_pos == k[j])[0][0]
                    min_E_pos[j] = min(min_E_pos[j], f_pos[indx])
                    max_E_pos[j] = max(max_E_pos[j], f_pos[indx])
                else:
                    indx = np.where(x_neg == k[j])[0][0]
                    min_E_neg[j] = min(min_E_neg[j], abs(f_neg[indx]))
                    max_E_neg[j] = max(max_E_neg[j], abs(f_neg[indx]))
        else:
            if abs_b: E = abs(E)
            if plot: plt.plot(k, E[i,:])
            min_E = np.minimum(E[i,:], min_E)
            max_E = np.maximum(E[i,:], max_E)
    # averaged spectrum
    stat_E = np.trapz(E[indt:,:], t[indt:], axis=0)/(t[-1] - t[indt])
    if hel:
        min_E_pos[np.where(min_E_pos == 1e30)] = \
                abs(stat_E[np.where(min_E_pos == 1e30)])
        min_E_neg[np.where(min_E_neg == 1e30)] = \
                abs(stat_E[np.where(min_E_neg == 1e30)])
        max_E_pos[np.where(max_E_pos == 0)] = \
                abs(stat_E[np.where(max_E_pos == 0)])
        max_E_neg[np.where(max_E_neg == 0)] = \
                abs(stat_E[np.where(max_E_neg == 0)])
    else:
        min_E[np.where(min_E == 1e30)] = \
                abs(stat_E[np.where(min_E == 1e30)])
        if abs_b:
            max_E[np.where(max_E == 0)] = \
                    abs(stat_E[np.where(max_E == 0)])
        else:
            max_E[np.where(max_E == 0)] = \
                    (stat_E[np.where(max_E == 0)])
    if plot:
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('$k$')

    if hel: return min_E_pos, min_E_neg, max_E_pos, max_E_neg, stat_E
    else: return min_E, max_E, stat_E

def local_max(k, E, order=1):

    """
    Function that computes the local maxima of a spectrum.

    Arguments:
        k -- array of wave numbers
        E -- spectrum E
        order -- order of the local maximum solver, which uses
                 scipy.signal.argrelextrema

    Returns:
        kmax -- position of the local maxima
        Emax -- values of the local maxima
    """

    from scipy.signal import argrelextrema

    inds_model_max = argrelextrema(E,
                 np.greater, order=order)
    kmax = k[inds_model_max]
    Emax = E[inds_model_max]

    return kmax, Emax

def compute_yks(k, E, N):

    """
    Function that computes N power law fittings, logarithmically
    equidistributed in k, of the spectrum E.

    Arguments:
        k -- array of wave numbers
        E -- array of spectrum values
        N -- number of power law fittings to discretize E
    """

    # compute N number of single power law fits around the model
    kps = np.logspace(np.log10(k[0]),
                      np.log10(k[-1]), N + 1)
    Ess = np.interp(kps, k, E)
    akss = np.zeros(N)
    c = np.zeros(N)
    akss[0], c[0] = slope(Ess[0], Ess[1],
                       kps[0], kps[1])
    kss = np.logspace(np.log10(kps[0]),
                      np.log10(kps[1]), 5)
    Ekss = kss**akss[0]*10**c[0]
    for i in range(2, N + 1):
        akss[i - 1], c[i - 1] = slope(Ess[i - 1], Ess[i],
                               kps[i - 1], kps[i])
        ksss = np.logspace(np.log10(kps[i - 1]),
                           np.log10(kps[i]), 5)
        Eksss = ksss**akss[i - 1]*10**c[i - 1]
        kss = np.append(kss, ksss)
        Ekss = np.append(Ekss, Eksss)
    #km, Em = mean_pos_loglog(np.append(kps[0], kps[2:]),
    #                         np.append(Ess[0], Ess[2:]))
    km, Em = mean_pos_loglog(kps, Ekss)

    return kss, Ekss, akss, km, Em, kps, c

def mean_pos_loglog(k, E):

    """
    Function that computes the loglog middle values km, EM of the intervals
    of the arrays k, E

    Arguments:
        k -- array of wave numbers
        E -- array of spectrum values

    Returns:
        km -- array of middle log values of the k intervals
        Em -- array of middle log values of the E intervals
    """

    N = len(k)
    km = np.zeros(N - 1)
    Em = np.zeros(N - 1)
    for i in range(1, N):
        km[i - 1] = np.sqrt(k[i - 1]*k[i])
        Em[i - 1] = np.sqrt(E[i - 1]*E[i])

    return km, Em

def slope(y1, y2, x1, x2):

    """
    Function that computes the slope between points 1 and 2

    Arguments:
        x1 -- x coordinate of point 1
        x2 -- x coordinate of point 2
        y1 -- y coordinate of point 1
        y2 -- y coordinate of point 2

    Returns:
        a -- slope between points 1 and 2
        c -- y-intercept of the straight line joining points 1 and 2
    """

    a = np.log10(y1/y2)/np.log10(x1/x2)
    c = np.log10(y1) - a*np.log10(x1)

    return a, c

def slope_A(x1, y1, x2, y2):

    """
    Function that computes the slope between points 1 and 2

    Arguments:
        x1 -- x coordinate of point 1
        x2 -- x coordinate of point 2
        y1 -- y coordinate of point 1
        y2 -- y coordinate of point 2

    Returns:
        a -- slope between points 1 and 2
        A -- amplitude of the fit y = A x^a
    """

    # slope
    a = np.log10(y1/y2)/np.log10(x1/x2)
    # amplitude
    A = y2*x2**(-a)

    return a, A

def red_blue_func(x, f, col=0):

    """
    Function that splits an array into positive and negative values, and
    assigns colors (red to positive and blue to negative).

    Arguments:
        x -- array of x
        f -- array of the function values
        col -- option to choose blue and red (default 0 is red for positive
               and blue for negative, 1 is swapped)

    Returns:
        x_pos -- array of x values where f is positive
        x_neg -- array of x values where f is negative
        f_pos -- array of f values where f is positive
        f_neg -- array of f values where f is negative
        color -- array of colors assigned (blue and red)
    """

    N = len(f)
    color = []
    f_pos=[]; x_pos=[]
    f_neg=[]; x_neg=[]
    for i in range(0, N):
        sgn = np.sign(f[i])
        if sgn > 0:
            if col == 0: color.append('red')
            if col == 1: color.append('blue')
            f_pos.append(f[i])
            x_pos.append(x[i])
        else:
            if col == 0: color.append('blue')
            if col == 1: color.append('red')
            f_neg.append(f[i])
            x_neg.append(x[i])
    f_pos = np.array(f_pos)
    f_neg = np.array(f_neg)
    x_pos = np.array(x_pos)
    x_neg = np.array(x_neg)

    return x_pos, x_neg, f_pos, f_neg, color

def plot_neg_pos(x, f, ls1='solid', lw1=1, ls2=':', lw2=2, col='black'):

    """
    Function that splits an array into positive and negative values, and
    plots them with different line styles.

    Arguments:
        x -- array of x
        f -- array of the function values
        col -- option to choose blue and red (default 0 is red for positive
               and blue for negative, 1 is swapped)
    """

    import matplotlib.pyplot as plt

    # plot positive and negative values with different line styles
    sgn = np.sign(f)
    converge = False
    sgn0 = sgn[0]
    i = 0
    lw = 1
    while not converge:
        sign = False
        i0 = i
        while not sign and not converge:
            if sgn0 == 1:
                ls = ls1
                lw = lw1
            else:
                ls = ls2
                lw = lw2
            if i==len(sgn) - 2: converge=True
            if sgn[i] != sgn0:
                sign = True
                sgn0 = sgn[i]
            i += 1
        plt.plot(x[i0-1:i], abs(f[i0-1:i]),
                 color=col, ls=ls, lw=lw)

def str_exp(exp, ak, den, diff=0.05):

    """
    Function that returns a string k^(a/den) if the absolute difference between
    the value a/den and the exponent ak is below diff.

    Arguments:
        exp -- initial string (given by the previous best estimation)
        ak -- slope
        den -- denominator of the fractions to be tested
        diff -- difference used to accept the corrected fraction approximating
                the slope

    Returns:
        exp -- updated string of the fractional slope
        diff -- difference updated
    """

    test = np.array(range(1, int(20*den)))
    test_i = test/den
    difft = abs(test_i - abs(ak))
    ind_min = np.argmin(difft)
    if difft[ind_min] < diff:
        m = test[ind_min]
        if ak > 0:
            if den == 1: exp = '$k^{%i}$'%m
            else: exp = '$k^{%i/%i}$'%(m, den)
        else:
            if den == 1: exp = '$k^{-%i}$'%m
            else: exp = '$k^{-%i/%i}$'%(m, den)
        diff = difft[ind_min]

    return exp, diff

def combine(k1, k2, E1, E2, facM, klim=10, exp=2):

    """
    Function that combines the spectra and wave number of two runs and uses
    the ratio between their magnetic amplitudes (facM) to compensate the
    GW spectrum by facM^2.

    Arguments:
        k1, k2 -- wave number arrays of runs 1 and 2
        E1, E2 -- GW spectral values arrays of runs 1 and 2
        facM -- ratio of the magnetic spectra amplitudes A2/A1
        klim -- wave number at which we switch from run2 to run 1
                (default 10)
        exp -- exponent used in facM to compensate the spectra (default 2,
               which correspond to that for GW spectra compensated by ratio
               between magnetic spectra)

    Returns:
        k -- combined wave number array
        E -- combined spectra
    """

    k = np.append(k2[np.where(k2 <= klim)], k1[np.where(k1 > klim)])
    E = np.append(E2[np.where(k2 <= klim)]/facM**exp, E1[np.where(k1>klim)])

    return k, E

def slopes_loglog(k, E):

    """
    Function that computes numerically the power law slope of a function
    E(k), taken to be the exponent of the tangent power law, i.e.,
    (\partial \ln E)/(\partial \ln k)

    Arguments:
        k -- independent variable
        E -- dependent variable

    Returns:
        slopes -- slope of E at each k in a loglog plot

    """
    slopes = np.zeros(len(k))
    slopes[0] = (np.log10(E[1]) - np.log10(E[0]))/ \
                        (np.log10(k[1]) - np.log10(k[0]))
    slopes[1] = (np.log10(E[2]) - np.log10(E[0]))/ \
                        (np.log10(k[2]) - np.log10(k[0]))
    for i in range(2, len(k) - 2):
         slopes[i] = (np.log10(E[i + 2]) + np.log10(E[i + 1]) \
                            - np.log10(E[i - 2]) - np.log10(E[i - 1]))/\
                            (np.log10(k[i + 1])+ \
                            np.log10(k[i + 2])-np.log10(k[i - 1]) - \
                            np.log10(k[i - 2]))
    slopes[-1] = (np.log10(E[-1]) - np.log10(E[-2]))/\
                        (np.log10(k[-1]) - np.log10(k[-2]))
    slopes[-2] = (np.log10(E[-1]) - np.log10(E[-3]))/\
                        (np.log10(k[-1]) - np.log10(k[-3]))
    return slopes

def get_min_max(f, E_a, E_b):

    """
    Function that returns the minimum and maximum of the power laws constructed
    for a different range of slopes.

    Arguments:
        f -- array of frequencies
        E_a -- 2d array of the minimum amplitudes (first index correspond to
               the slope and the second to the frequency)
        E_b -- 2d array of the power laws corresponding to the maximum
               amplitudes

    Returns:
        minE -- array of minimum values of the spectra over all slopes
        maxE --  array of maximum values of the spectra over all slopes
    """

    minE = np.zeros(len(f)) + 1e30
    maxE = np.zeros(len(f))
    for i in range(0, len(f)):
        good = np.where(E_a[:, i] != 0)
        minE[i] = np.min(E_a[good, i])
        maxE[i] = np.max(E_b[:, i])

    return minE, maxE
