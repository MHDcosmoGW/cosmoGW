"""
run.py is a Python routine that contains the class run used to store all the
variables associated to a specific run, as well as functions to initialize
and postprocess a set of runs given by an array of directories.

dirs.py contains lists of directories for specific sets of runs in the
literature.

The runs correspond to Pencil Code simulations of gravitational waves from
different sources in the early universe, for example, from MHD turbulence,
see cosmoGW project (https://github.com/AlbertoRoper/cosmoGW)
for details.

Author: Alberto Roper Pol
Created: 01/01/2021
Updated: 01/11/2023 (release of the cosmoGW code)
"""

import os
HOME = os.getcwd() + '/'
import reading as re
import warnings
from dirs import read_dirs as rd
import numpy as np
import pandas as pd
import pickle
import spectra as spec

def initialize_runs(R=[], dir0=HOME, dirs='', quiet=True, opt=0, debug=False, project='',
                    char=False, load=False, ext=''):

    """
    Function to initialize the dictionary with the list of runs pointing to
    the variables of the class run, which contains the spectra, the time
    series, and direct calculations.

    Arguments:
        R -- array with the name of the runs to be read
        dir0 -- directory that contains the runs to be read
        dirs -- dictionary with the name of each of the directories to be read
        quiet -- prints the list of read runs if False (default True)
        opt -- option to choose some reading routines (default is 0 and 1 can
               be chosen if 0 gives warnings)
        debug -- option to print out debugging tests (default is False)
        project -- name of the specific project (if dirs, R are not given, they can be
                   read from the specific project if it exists in dirs.py)
        char -- option to execute postprocessing for each run using the function
                'characterize_runs' (default False)
        load -- option to load the runs from previously generated pickle files, which
                already contains the results from the simulations (faster than reading
                from the files, but this needs to be done at least once)
        ext -- possible extension to select a subset of runs of the specified project

    Returns:
        runs -- dictionary with the initialized values of the runs
    """

    # if dirs is not given, it looks for it in dirs.py from the specific project
    if dirs=='':
        if  project != '':
            dir0 = HOME + project + '/'
            dirs = {}
            dirs = rd(project, dirs, ext=ext)
            R = [s for s in dirs]
            print("\n dirs is being set up for project: ", project)
            if ext != '': print('extension: ', ext)
            print('\n', dirs, '\n')
        else:
            print('A valid project name needs to be given to initialize_runs function!')
    else:
        print('dirs is read from the input (needs to be correct)')
        print(dirs)

    runs = {}
    for i in R:

        dir_run = dirs.get(i)
        rd_pick = True

        if not load:
            print('Reading files from the simulations')
            run_var = run(i, dir0, dir_run, quiet=quiet, opt=opt,
                          debug=debug, project=project)

        else:
            print('Reading pickle variables generated from simulations: ', i)
            try:
                f = open(dir0 + dir_run + '/' + i + '.pckl', 'rb')
                run_var = pickle.load(f)
                f.close()
            except:
                rd_pick = False
                print('reading pickle file for ', i, ' failed! make sure you first initialize and',
                      ' generate the pickle files from the run directories')
        if rd_pick:
            if char: run_var.characterize_run(quiet=quiet)
            runs.update({i:run_var})

    if not quiet:
        print('\n The runs that have been read are:')
        print([s for s in runs.keys()])

    return runs

def characterize_runs(runs, quiet=True):

    """
    Function that executes the characterize_run function contained
    within the class 'run' for each run in the dictionary 'runs' that
    contains each of the runs to be analyzed.

    Arguments:
        runs -- dictionary of variables of the class run
        quiet -- prints the variables if False (default True)

    Returns:
        runs -- updated dictionary of variables of the class run
    """

    for i in runs:
        run = runs.get(i)
        run.characterize_run(quiet=quiet)
        
def save_runs(runs, dir0=''):

    """
    Function that saves the run variables as pickle variables.

    Arguments:
        runs -- dictionary of the run variables
        dir0 -- directory where the pickle variables are saved (by default
                same as the run directory)
    """

    for i in runs:
        run = runs.get(i)
        run.save(dir0=dir0)

class run():

    """
    Class that contains the variables corresponding to a single run of the
    Pencil Code of GW generation from early universe sources.
    """

    def __init__(self, name_run, dir0, dir_run, quiet_war=False, quiet=True,
                 opt=0, debug=False, project=''):

        """
        Initialize class run and reads the spectra and the time series
        from the files power_#sp.dat and time_series.dat in the directory
        dir0 + dir_run + '/data/', where #sp corresponds to a specific
        spectrum.

        For a list of the spectra #sp, run: print(run.spectra_avail)

        Arguments:
            name_run -- name used to identify the specific run
            dir0 -- directory that contains all the runs to be read
            dir_run -- specific directory of the run
            quiet_war -- used to ignore warnings if set to True (default False)
            quiet -- prints the list of the run spectra if False (default True)

        Returns:
            run -- initialized class run with the variables:
                   name_run -- name used to identify the specific run
                   dir_run -- specific directory of the run
                   spectra -- directory that includes the different spectra
                   spectra_avail -- list of available spectra within the
                                    spectra dictionary
                   ts -- dictionary that includes the different variables
                         in the time series data
                   ts_avail -- list of available variables within the ts
                               dictionary
                   df -- data frame with info of the run
                   proj -- name of the specific project that the runs belong to
                           (usually same as dir_run)
        """

        # Option to ignore warnings
        if quiet_war:
            warnings.filterwarnings('ignore',
                                       category=np.VisibleDeprecationWarning)
        else:
            warnings.filterwarnings('error',
                                       category=np.VisibleDeprecationWarning)

        self.name_run = name_run
        self.dir_run = dir_run
        print('Reading run ' + name_run + '\n')

        # Reading spectra
        print(dir0, self.dir_run)
        self.spectra = re.read_spectra_runs(dir_data=dir0 + self.dir_run + '/data/', opt=opt)
        keys = self.spectra.keys()
        self.spectra_avail = [s for s in self.spectra.keys() if not s=="k"]
        self.spectra_avail = [s for s in self.spectra_avail if not s=="k0"]
        self.spectra_avail = [s for s in self.spectra_avail if not 't_' in s]
        if not quiet:
            print('Spectra computed: ', self.spectra_avail)
            print('\n')

        # Reading time series
        self.ts = re.read_ts(dir_data=dir0 + self.dir_run + '/data/',
                             opt=opt, debug=debug)
        self.ts_avail = [s for s in self.ts.keys() if not s=='it']
        self.ts_avail = [s for s in self.ts_avail  if not s=='t']
        self.ts_avail = [s for s in self.ts_avail  if not s=='dt']
        if not quiet:
            print('Time series computed: ', self.ts_avail)
            print('\n')
        
        self.proj = project
        self.df = pd.DataFrame({'project': [self.proj], 'dir': [self.dir_run],
                                'name': [self.name_run]})

    def characterize_run(self, quiet=True, upd_Pi=True, upd_EGW=True, comp_EGW_tot=True,
                         comp_integ_max=True, comp_max_fields=True, ts_max=['EEM', 'EEK'],
                         sp_max=['mag', 'kin'], max_allowed_Om=1, cs2=1/3,
                         comp_rho=False, comp_max_en=False, fields_EE=['EEM', 'EEK'],
                         max_fields=['b', 'u']):

        """
        Function that computes the results used to characterize the run
        using the spectra and the time series.

        Arguments:
            quiet -- option to avoid printing some info of the runs (default True)
                     using print_characterize
            upd_Pi -- option to compute the Pi spectrum using update_Pi (default True)
            upd_EGW -- option to compute the GW spectrum EGW and OmGW using update_EGW
                       (default True)
            comp_EGW_tot -- option to combine fields to compute the actual EGW and OmGW
                            using update_EGW (default True)
            comp_integ_max -- option to compute integrated values of each of the spectra,
                              as well as the maximum Emax and kpeak as a function of time
                              (default True)
            comp_max_fields -- option to compute the maximum values of averaged fields
                               used to characterize the type of simulation (default True)
            ts_max -- time series of averaged fields computed using comp_max_fields
                      (default is 'EEM' and 'EEK', i.e., magnetic and velocity fields, used
                      for MHD simulations)
            sp_max -- spectra of fields used to compute the maximum values of the fields
                      (default is 'mag' and 'kin', i.e., magnetic and velocity spectra, used
                      for MHD simulations)
            max_allowed_Om -- maximum value allowed of averaged fields to show a warning
                              (default is 1 for normalized fields in the early universe)
            cs2 -- speed of sound (eos) used to compute the characteristic speeds (e.g.,
                   Alfvenic speed for magnetic fields) (default is 1/3, i.e., RD era)
            comp_rho -- option to compute the energy density (default True)
            comp_max_en -- option to compute the max energies (default True)
            fields_EE -- values of the time series where the energies are stored to compute
                         the maximum energy (if comp_max_en is True) (default are 'EEM' and 'EEK'
                         for MHD turbulence)
            max_fields -- values of the time series containing the maximum over the domain
                          (if comp_max is True) (default are 'b' and 'u' for MHD turbulence)
        """

        # compute the spectrum of the stresses Pi = Str/k^2 and the GW spectra
        # these functions are adapted for RD runs (needs to be adapted for
        # other cases)
        if upd_Pi: self.update_Pi()
        if upd_EGW: self.update_EGW(comp_tot=comp_EGW_tot)

        # compute maximum and integrated values of spectra and spectral peaks
        if comp_integ_max: self.compute_integ_max_spectra()

        # compute the time at which the characteristic fields (mag and kin by default)
        # have their maxima, the max value and the position of the spectral peak at
        # that time
        if comp_max_fields:
            Nts_max = len(ts_max)
            if len(sp_max) != Nts_max:
                print('provide same number of ts_max and sp_max in function characterize_run')
            else:
                self.Om_max = np.zeros(Nts_max)
                self.t_max = np.zeros(Nts_max)
                self.kf_max = np.zeros(Nts_max)
                self.v_char = np.zeros(Nts_max)
                self.t_eddy = np.zeros(Nts_max)
                for i in range(0, Nts_max):
                    self.Om_max[i], self.t_max[i], self.kf_max[i] = \
                        self.check_max_spectra_ts(sp_max[i], ts_max[i])
                    # set maximum allowed Om to consider it a correct value and show a
                    # warning if the value is over the allowed one
                    if self.Om_max[i] > max_allowed_Om:
                         print('Maximum value of the ', ts_max[i], ' energy density is too',
                              ' large: ', self.OmMmax[i], ' > ', max_allowed_Om)
                    # characteristic speed based on total enthalpy
                    # (e.g., Alfven speed for magnetic field) and turnover time
                    self.v_char[i] = np.sqrt(2*self.Om_max[i]/(1 + cs2))
                    if self.v_char[i] != 0 and self.kf_max[i] != 0:
                        self.t_eddy[i] = 2*np.pi/self.kf_max[i]/self.v_char[i]

                # nature of the source is given by the field that has its maximum
                ind_max = np.argmax(self.Om_max)
                self.type = ts_max[ind_max]
                self.Ommax = self.Om_max[ind_max]

        # compute energy density time series
        if comp_rho: self.compute_rho()
        
        # compute total turbulent and maximum energy densities from time
        # series
        if comp_max_en: self.compute_total_max_energies(ts_fields=fields_EE,
                                                        max_fields=max_fields)
        
        # print some info characterizing the run when quiet is False
        if not quiet and comp_max_fields:
            print(self.name_run, '(', self.type, '): Omega max: ', self.Ommax,
                  ', kf: ', self.kf_max[ind_max], 'v: ', self.v_char[ind_max],
                  'te: ', self.t_eddy[ind_max])

    def update_Pi(self, case_run='RD'):

        """
        Function that computes the spectrum Pi (k, t) as the spectrum of the
        projected stress (Str) divided by k^2 (same for helical Str).
        
        Arguments:
            case_run -- specific type of run (used to compensate Str), see below
                        (default is RD)

        The spectra dictionary of run is updated with:
            Pi -- spectrum Str divided by k^2
            t_Pi -- time array of spectrum Pi
            helPi -- spectrum helStr divided by k^2
            t_helPi -- time array of spectrum helPi
            
        Note that the Str is computed from the source of the GW equation in the Pencil
        Code, so it includes an additional 6/t factor when dealing with RD era computations.
        
        We then need to multiply the Str spectrum by t^2/36 to compensate for this prefactor
        and get the physical spectrum. At the moment, this is corrected only for RD calculations,
        to be updated for other cases!
        
        An additional factor of 2 is included since we define Str from T+^2 + Tx^2 in the Pencil Code,
        such that Tij Tij = 2 (T+^2 + Tx^2), used to define Pi
        
        Reference: A. Roper Pol, A. Brandenburg, T. Kahniashvili, A. Kosowsky, S. Mandal,
        "The timestep constraint in solving the gravitational wave equations sourced by hydromagnetic
        turbulence," Geophys. Astrophys. Fluid Dynamics 114, 1, 130 (2020), arXiv:1807.05479,
        eq. A8 for 6/t factor, eq. B14 for general definition of GW and Str spectra.
        """

        # check if Str is within the available spectra
        if 'Str' in self.spectra_avail:

            k = self.spectra.get('k')
            t = self.spectra.get('t_Str')
            sp = self.spectra.get('Str')
            tij, kij = np.meshgrid(t, k, indexing='ij')
            # avoid division by 0
            kij[np.where(kij == 0)] = 1.

            if case_run == 'RD': pref = 2*tij**2/36
            else: print('no other options besides RD are implemented to compute Pi')

            Pi = sp/kij**2*pref
            Pi[np.where(kij == 0)] = 0.
            self.spectra.update({'Pi': Pi})
            self.spectra.update({'t_Pi': t})
            
            if 'Pi' not in self.spectra_avail:
                self.spectra_avail.append('Pi')

            # check if helStr is within the available spectra
            if 'helStr' in self.spectra_avail:

                t = self.spectra.get('t_helStr')
                tij, kij = np.meshgrid(t, k, indexing='ij')

                if case_run == 'RD': pref = 2*tij**2/36
                else: print('no other options besides RD are implemented to compute Pi')

                kij[np.where(kij == 0)] = 1.
                sphel = self.spectra.get('helStr')
                helPi = sphel/kij**2*pref
                helPi[np.where(kij == 0)] = 0.
                self.spectra.update({'helPi':helPi})
                self.spectra.update({'t_helPi':t})

                if 'helPi' not in self.spectra_avail:
                    self.spectra_avail.append('helPi')
                    
    def update_EGW(self, case_run='RD', comp_tot=True):

        """
        Function that computes the spectra of the GW energy density in two ways:
        1) only from the spectra of the time derivative of the strains (incomplete
        description in the early universe) but correct for later times and hence, for
        observational prospects, 2) including mixed terms and strains spectra to
        take into account the terms coming from the difference between time derivatives
        with respect to cosmic and conformal time.
        
        Reference: A. Roper Pol, A. Brandenburg, T. Kahniashvili, A. Kosowsky, S. Mandal,
        "The timestep constraint in solving the gravitational wave equations sourced by hydromagnetic
        turbulence," Geophys. Astrophys. Fluid Dynamics 114, 1, 130 (2020), arXiv:1807.05479,
        appendix B
        
        An additional factor of 2 is included since we define GW spectra from h+^2 + hx^2 in the Pencil Code,
        such that hij hij = 2 (h+^2 + hx^2), used to define, for example, GWh.
        This factor, times the 1/12 obtained from the normalization valid during RD era gives a prefactor of 6.
        
        At the moment, this is corrected only for RD calculations,to be updated for other cases!

        The spectra dictionary of run are updated with:

            EGW -- GW energy density spectrum (linear interval) using only GWs
            OmGW -- GW energy density spectrum (logarithmic interval)
            EGW_tot -- GW energy density spectrum (linear interval)
                       using GWs and GWh and/or GWm
            OmGW_tot -- GW energy density spectrum (logarithmic interval)
                       using GWs and GWh and/or GWm
            helEGW -- helical analogous spectrum to EGW
            helOmGW -- helical analogous spectrum to OmGW
            helEGW_tot -- helical analogous spectrum to EGW_tot
            helOmGW_tot -- helical analogous spectrum to OmGW_tot
            t_#sp -- corresponding time arrays of each #sp
        """
        
        if case_run == 'RD': pref = 1/6
        else: print('no other options besides RD are implemented to compute Pi')

        k = self.spectra.get('k')
        
        # function used below to read the spectrum of the time derivative strains
        # and compute EGW = pref*GWs and OmGW = EGW*k
        
        def read_GW(k, sp_n='GWs', hel=''):
            
            t = self.spectra.get('t_' + hel + sp_n)
            tij, kij = np.meshgrid(t, k, indexing='ij')
            sp = self.spectra.get(hel + sp_n)*pref
            
            self.spectra.update({hel + 'EGW': sp})
            self.spectra.update({hel + 'OmGW': kij*sp})
            self.spectra.update({'t_' + hel + 'EGW': t})
            self.spectra.update({'t_' + hel + 'OmGW': t})
            
            if (hel + 'EGW') not in self.spectra_avail:
                self.spectra_avail.append(hel + 'EGW')
            if (hel + 'OmGW') not in self.spectra_avail:
                self.spectra_avail.append(hel + 'OmGW')
                
            return kij, t, tij, sp
            
        # function used below to read the spectrum of the time derivative strains (GWs),
        # the strains (GWh), and the mixed term (GWm), and compute
        # EGW_tot = pref*(GWs + GWh/t^2 - 2 GWm/t) and OmGW_tot = EGW_tot*k

        def read_GW_comb(kij, t, tij, sp1, sp2_n='GWh', sp3_n='GWm', hel=''):
            
            extra = False
            if (hel + sp2_n) in self.spectra_avail:
                extra = True
                sp2 = self.spectra.get(hel + sp2_n)

                # check if shape (Nk, Nt) of spectra is the same as the one sp1
                # and if not, cuts spectra GWh and GWm
                sh_sp2 = np.shape(sp2)[0]
                if np.shape(sp1)[0] > sh_sp2:
                    print(hel + 'GWs has more time steps than' + hel + 'GWh, check!')
                    sp1 = sp1[:sh_sp2, :]
                    tij = tij[:sh_sp2, :]
                    kij = kij[:sh_sp2, :]
                    t = t[:sh_sp2]
                sp1 += sp2/tij**2*pref
            
            if (hel + sp3_n) in self.spectra_avail:
                extra = True
                sp3 = self.spectra.get(hel + sp3_n)
                sp1 -= 2*sp3/tij*pref
                
            if extra:
                self.spectra.update({hel + 'EGW_tot': sp1})
                self.spectra.update({hel + 'OmGW_tot': kij*sp1})
                self.spectra.update({'t_' + hel + 'EGW_tot': t})
                self.spectra.update({'t_' + hel + 'OmGW_tot': t})

                if (hel + 'EGW_tot') not in self.spectra_avail:
                    self.spectra_avail.append(hel + 'EGW_tot')
                if (hel + 'OmGW_tot') not in self.spectra_avail:
                    self.spectra_avail.append(hel + 'OmGW_tot')
        
        # read spectra of time derivatives of strains h' (GWs)
        if 'GWs' in self.spectra_avail:
            kij, t, tij, sp = read_GW(k, sp_n='GWs', hel='')
            if comp_tot:
                # if GWh (spectra of strains h) and/or GWm (mixed spectra of h and h')
                # are computed, then combine them to compute the total EGW during the simulation
                read_GW_comb(kij, t, tij, sp, sp2_n='GWh', sp3_n='GWm', hel='')

        # same calculation for chiral part of the GW spectrum
        if 'helGWs' in self.spectra_avail:
            kij, t, tij, sp = read_GW(k, sp_n='GWs', hel='hel')
            if comp_tot:
                read_GW_comb(kij, t, tij, sp, sp2_n='GWh', sp3_n='GWm', hel='hel')

    def compute_integ_max_spectra(self):

        """
        Function that computes the integrated and maximum values of the spectra, as
        well as the position of the spectral peak.

        Updates the content of the variable run with:
            #sp_max   -- time dependent maximum value of each spectral function
            #sp_integ -- time dependent integrated value of the spectrum
            #sp_kpeak -- time dependent spectral peak
        """

        k = self.spectra.get('k')
        for m in self.spectra_avail:
            t = self.spectra.get('t_' + m)
            Nt = len(t)
            kpeak = np.zeros(Nt)
            Emax = np.zeros(Nt)
            Emean = np.zeros(Nt)
            E = self.spectra.get(m)
            for i in range(0, Nt):
                kpeak[i], Emax[i] = spec.compute_kpeak(k[1:], E[i, 1:], quiet=True)
                Emean[i] = np.trapz(E[i, :], k)
            self.spectra.update({m + '_kpeak': kpeak})
            self.spectra.update({m + '_max': Emax})
            self.spectra.update({m + '_integ': Emean})

    def check_max_spectra_ts(self, E, EE):

        """
        Function that computes the maximum value of a field over time.
        It computes the maximum of the integrated spectrum E and the maximum
        of the averaged field EE (from the time series) and returns the largest
        value of the two.

        Arguments:
            E -- string that defines the name of the spectral function
                 (e.g., E = 'mag')
            EE -- string that defines the name of the time series variable
                  (e.g., EE = 'EEM')

        Returns:
            max -- maximum over time of the averaged value
            tmax -- time at which the field mean value is maximum
            kf -- spectral peak at tmax
        """

        EEmmax = 0
        t = self.ts.get('t')
        tmax = t[0]
        kf = 0

        if EE in self.ts_avail:

            EEm = self.ts.get(EE)
            indmax_ts = np.argmax(EEm)
            EEmmax = EEm[indmax_ts]
            tmax = t[indmax_ts]

        if E in self.spectra_avail:

            # compute integrated spectrum if it has not been computed
            if (E + '_integ') not in self.spectra:
                self.compute_integ_max_spectra()

            mean = self.spectra.get(E + '_integ')
            t = self.spectra.get('t_' + E)
            indmax = np.argmax(mean)
            kpeak = self.spectra.get(E + '_kpeak')
            if mean[indmax] > EEmmax:
                EEmmax = mean[indmax]
                tmax = t[indmax]
                kf = kpeak[indmax]
            else:
                if EE in self.ts_avail:
                    kf = np.interp(self.ts.get('t')[indmax_ts], t, kpeak)

        return EEmmax, tmax, kf
                    
    def compute_rho(self):

        """
        Function that computes the energy density if urms and EEK are available
        in the time series.

        The ts dictionary of run is updated with:
            rho -- energy density
        """

        if 'rho' not in self.ts_avail:
            if 'urms' in self.ts_avail and 'EEK' in self.ts_avail:
                urms = self.ts.get('urms')
                good = np.where(urms != 0)
                EEK = self.ts.get('EEK')
                rho = EEK**0
                rho[good] = 2*EEK[good]/urms[good]**2
                self.ts.update({'rho': rho})
                self.ts_avail.append('rho')

    def compute_total_max_energies(self, ts_fields=['EEM', 'EEK'], max_fields=['u', 'b']):

        """
        Function that updates run with the time series of the total turbulent
        energy density 'EEtot', and maximum kinetic 'EEKmax',
        magnetic 'EEMmax', and total 'EEtotmax' energy densities,
        when they are computed.

        The ts dictionary of run is updated with:
            EEtot -- total turbulent energy density (EEM + EEK)
            EEKmax -- maximum kinetic energy density
            EEMmax -- maximum magnetic energy density
            EEtotmax -- sum of maximum kinetic and magnetic energies
        """

        # add up over the fields described in ts_fields to compute the maximum energy
        # (e.g. kinetic + magnetic energy density)
        if 'EEtot' not in self.ts_avail:        
            j = 0
            for i in ts_fields:
                if i in self.ts_avail:
                    EE = self.ts.get(i)
                    if j == 0: tot_E = EE
                    else: tot_E += EE
                    j += 1
            # if some field is available, we return total energy density
            if j > 0:
                self.ts.update({'EEtot': tot_E})
        
        # compute max energy density from fields in max_fields if available
        for i in max_fields:
            if ('EEmax_' + i) not in self.ts_avail:
                if (i + 'max') in self.ts_avail:
                    EEmax = self.ts.get(i + 'max')**2/2
                    self.ts.update({'EEmax_' + i: EEmax})
                    self.ts_avail.append('EEmax' + i)

    def compute_pol(self, sp='EGW', A=1., exp=0):

        """
        Function that computes the GW and magnetic polarization spectra
        PGW(k), Ph(k), and PM(k).
        
        Arguments:
            sp -- spectra used to compute the spectrum of polarization
                  (antisymmetric over symmetric spectra)
            A, exp -- values to compensate the antisymmetric spectrum with
                      A k^exp (e.g., A = .5, exp = 1 for magnetic helicity)
        """
        
        sp_hel = 'hel' + sp
        P_sp = 'P_' + sp
        if P_sp not in self.spectra_avail:
            if sp in self.spectra_avail and sp_hel in self.spectra_avail:
                EE = self.spectra.get(sp)
                hel_EE = self.spectra.get(sp_hel)
                t = self.spectra.get('t_' + sp_hel)
                t2 = self.spectra.get('t_' + sp)
                if len(t) != len(t2):
                    print('spectrum ', sp, ' and ', sp_hel, ' of run ', self.name_run, 'do not have same number of',
                          ' data points in time')
                    inds = range(0, min(len(t), len(t2)))
                    EE = EE[inds, :]
                    hel_EE = hel_EE[inds, :]
                    t = t[inds]
                Ak = A
                if exp != 0:
                    k = self.spectra.get('k')
                    tij, kij = np.meshgrid(t, k, indexing='ij')
                    Ak *= kij**exp
                    Ak[kij == 0] = 1e-20
                PP = np.zeros((np.shape(EE)))
                inds = np.where(EE != 0)[0]
                PP[inds] = Ak*hel_EE[inds]/EE[inds]
                self.spectra.update({P_sp: PP})
                self.spectra.update({'t_' + P_sp: t})
                self.spectra_avail.append(P_sp)

    def min_max_stat(self, abs_b=True, sp='GWs', indt=0, indtf=-1,
                     plot=False, hel=False):

        """
        Function that computes the minimum, the maximum, and the
        averaged functions over time of a spectral function of the run.

        Arguments:
            abs_b -- option to take absolute value of spectrum (default True)
            sp -- string that indicates the spectral function
            indt, indtf -- indices of time array to perform the average
                           from t[indt] to t[indtf] (default is 0 to final time)
            plot -- option to overplot all spectral functions
                    from t[indt] to t[indtf]
            hel -- option for helical spectral functions where positive and
                   negative values can appear (default False)
                   It then returns min_E_pos, min_E_neg, max_E_pos, max_E_neg
                   referring to the maximum/minimum absolute values of the positive
                   and negative values of the helical funtion.

        The updated spectral functions within the run.spectra dictionary
        are:
            #sp_min_sp -- minimum of the spectral function #sp over times
            #sp_max_sp -- maximum of the spectral function #sp over times
            #sp_stat_sp -- average of the spectral function #sp over times
                           from t[indt] to t[-1]
            if hel then also updates positive and negative values separately
                    #sp_pos_min, #sp_neg_min, #sp_pos_max, #sp_neg_max
        """

        if sp in self.spectra_avail:
            t = self.spectra.get('t_' + sp)
            E = self.spectra.get(sp)
            k = self.spectra.get('k')
            if hel:
                min_sp_pos, min_sp_neg, max_sp_pos, max_sp_neg, stat_sp = \
                        spec.min_max_stat(t, k, E, abs_b=abs_b, indt=indt,
                                          indtf=indtf, plot=plot, hel=True)
                self.spectra.update({sp + '_pos_min_sp':min_sp_pos})
                self.spectra.update({sp + '_neg_min_sp':min_sp_neg})
                self.spectra.update({sp + '_pos_max_sp':max_sp_pos})
                self.spectra.update({sp + '_neg_max_sp':max_sp_neg})
                self.spectra.update({sp + '_min_sp': \
                                     np.minimum(min_sp_pos, min_sp_neg)})
                self.spectra.update({sp + '_max_sp': \
                                     np.maximum(max_sp_pos, max_sp_neg)})
                self.spectra.update({sp + '_stat_sp':stat_sp})
            else:
                min_sp, max_sp, stat_sp = \
                        spec.min_max_stat(t, k, E, abs_b=abs_b, indt=indt,
                                             plot=plot)
                self.spectra.update({sp + '_min_sp': min_sp})
                self.spectra.update({sp + '_max_sp': max_sp})
                self.spectra.update({sp + '_stat_sp': stat_sp})

        else: print(sp, ' spectrum is not available!')

    def save(self, dir0=''):

        """
        Function that saves the run in a pickle variable.

        Arguments:
            dir0 -- directory where to save the variable
                    (default is the current directory)
        """

        # change to dir0 if it is given, otherwise stay current directory
        if dir0 == '': dirr = dir0 + self.dir_run
        else: dirr = dir0
        cwd = os.getcwd()
        os.chdir(dirr)

        name_f = self.name_run + '.pckl'
        f = open(name_f, 'wb')
        print('Saving ' + self.name_run + '\n')
        print('The output file is ' + name_f + '\n saved in the directory',
              os.getcwd() + '\n')
        pickle.dump(self, f)
        f.close()

        # return to initial directory
        os.chdir(cwd)
                
    ############################# NOT PUBLIC ON GITHUB ########################

    def compute_GWsmodel_Pi0(self, tfin=1e4):

        """
        Function that computes the GW energy density spectrum OmGW(k) using the
        initial time Pi spectrum, considering the assumption that Pi does not
        evolve in time.
        """

        import numpy as np
        import GW_analytical as an

        if 'Pi' in self.spectra_avail:
            k = self.spectra.get('k')
            t = self.spectra.get('t_Pi')
            #Pi = np.array(self.spectra.get('Pi')[:, 1:], dtype='float')
            Pi0 = self.spectra.get('Pi')[0, :]
            #for i in range(0, len(k)):
            #    Pi0[i] = np.interp(self.tini, t, Pi[:, i])
            EGW_num, OmGW_num = an.OmGW_from_Pi0(Pi0, k, t, tfin=tfin,
                                                 tini=self.tini)
            self.spectra.update({'EGW_num_Pi0': EGW_num})
            self.spectra.update({'OmGW_num_Pi0': OmGW_num})
            if not 'EGW_num_Pi0' in self.spectra_avail:
                self.spectra_avail.append('EGW_num_Pi0')
                self.spectra_avail.append('OmGW_num_Pi0')
            if 'PiM_model_EM' in self.spectra_avail:
                Pi_mod = self.spectra.get('PiM_model_EM')[0, :]
                ks_mod = self.spectra.get('k_PiM_model_EM')
                EGW_mod, OmGW_mod = an.OmGW_from_Pi0(Pi_mod, ks_mod, t,
                                                     tfin=tfin, tini=self.tini)
                self.spectra.update({'EGW_mod_Pi0': EGW_mod})
                self.spectra.update({'OmGW_mod_Pi0': OmGW_mod})
                if not 'EGW_mod_Pi0' in self.spectra_avail:
                    self.spectra_avail.append('EGW_mod_Pi0')
                    self.spectra_avail.append('OmGW_mod_Pi0')

    #def compute_PiM(self, ks=[], alpha=2, indt=0, alp_str='2'):

        #import GW_analytical as an
        #import numpy as np

        #if len(ks) == 0: ks = np.logspace(-2, 3, 10000)
        #if 'mag' in self.spectra_avail:
            # k = self.spectra.get('k')[1:]
            # mag = self.spectra.get('mag')[:, 1:]
            # # fit magnetic spectrum to smoothed broken power law at initial
            # # time
            # _ = an.model_EM(k, mag, ks=ks)
            # A_EM = _[3]
            # xc_EM = _[5]
            # # compute Pi from magnetic field parameters assuming Gaussianity
            # kf_Pi, PiM = an.compute_PiM(A_EM, xc_EM, alpha=alpha,
            #                             alp_str=alp_str)
            # PiM = 10**np.interp(np.log10(ks), np.log10(kf_Pi),
            #                     np.log10(PiM))
            # self.spectra.update({'Pi_model': PiM})
            # self.spectra.update({'ks_Pi_model': ks})
            # if not 'Pi_model' in self.spectra_avail:
            #     self.spectra_avail.append('Pi_model')

    # def compute_GWsmodel_coh(self, NNt=100):
    #
    #     import GW_analytical as an
    #     import numpy as np
    #
    #     if 'Pi' in self.spectra_avail:
    #         Pi = np.array(self.spectra.get('Pi'), dtype='float')
    #         k = self.spectra.get('k')
    #         t = self.spectra.get('t_Pi')
    #         EGW, OmGW = an.OmGW_from_Pi_coh(t, k, Pi, NNt=NNt)
    #         self.spectra.update({'EGW_model_coh': EGW})
    #         self.spectra.update({'OmGW_model_coh': OmGW})
    #         if not 'EGW_model_coh' in self.spectra_avail:
    #             self.spectra_avail.append('EGW_model_coh')
    #             self.spectra_avail.append('OmGW_model_coh')

    def compute_GWsmodel_coh_PiM(self, NNt=100):

        import GW_analytical as an
        import numpy as np

        if 'Pi' in self.spectra_avail:
            Pi = np.array(self.spectra.get('Pi'), dtype='float')
            k = self.spectra.get('k')
            t = self.spectra.get('t_Pi')
            EGW, OmGW = an.OmGW_from_Pi_coh(t, k, Pi,
                                            tini=self.tini,NNt=NNt)
            self.spectra.update({'EGW_num_coh': EGW})
            self.spectra.update({'OmGW_num_coh': OmGW})
            if not 'EGW_num_coh' in self.spectra_avail:
                self.spectra_avail.append('EGW_num_coh')
                self.spectra_avail.append('OmGW_num_coh')
            if 'PiM_model_EM' in self.spectra_avail:
                Pi_mod = np.array(self.spectra.get('PiM_model_EM'),
                                  dtype='float')
                ks_mod = self.spectra.get('k_PiM_model_EM')
                EGW_mod, OmGW_mod = an.OmGW_from_Pi_coh(t, ks_mod, Pi_mod,
                                                        tini=self.tini,NNt=NNt)
                self.spectra.update({'EGW_mod_coh': EGW_mod})
                self.spectra.update({'OmGW_mod_coh': OmGW_mod})
                if not 'EGW_mod_coh' in self.spectra_avail:
                    self.spectra_avail.append('EGW_mod_coh')
                    self.spectra_avail.append('OmGW_mod_coh')

    def check_initcond_mag_old(self, kpeak=1., a=4., b=5/3, alp=2, indt=0,
                               p0=1e-5, Nk=100000):

        import scipy.optimize as opt
        import GW_analytical as an
        import numpy as np

        # old fit for smoothed broken power law, ignoring scaling D1, D2 factors
        # that determines the initial spectrum for the vector potential
        # fix all parameters except the amplitude
        def fit_test(k, A):
            return an.fit_smoothed_bPL(k, A=A, a=a-2, b=b+2, alp=alp, xc=kpeak,
                                       lfactors=False)

        # fit amplitude using the original model for EM (k)
        # from the Pencil Code (note that it is defined for the vector
        # potential)
        k = self.spectra.get('k')
        if 'mag' in self.spectra_avail:
            tk = self.spectra.get('t_mag')
        elif 'GWs' in self.spectra_avail:
            tk = self.spectra.get('t_GWs')
        elif 'kin' in self.spectra_avail:
            tk = self.spectra.get('t_kin')
        self.df['tfin'] = tk[-1]
        self.df['t0'] = tk[0]
        self.df['k0'] = k[1]
        if 'mag' in self.spectra_avail:
            self.df['a'] = a
            self.df['b'] = b
            EM = self.spectra.get('mag')
            popt = opt.curve_fit(fit_test, k[1:], EM[indt, 1:]/k[1:]**2,
                                 p0=(p0))
            AA = popt[0][0]
            # use fit to generate high resolution spectrum
            ks = np.logspace(np.log10(k[1]), np.log10(k[-1]), Nk)
            Ak = fit_test(ks, A=AA)
            EM_mod = ks**2*Ak
            self.spectra.update({'EM0_model':EM_mod})
            self.spectra.update({'k_EM0_model':ks})

            # compute real values of the EM spectral peak and its position
            indkp = np.argmax(EM_mod)
            if indt == 0:
                self.df['EM0peak'] = EM_mod[indkp]
                self.df['kpeak'] = ks[indkp]

            # compute real alpha parameter using D1, D2 factors
            # (should be roughly the same!!) by fixing all parameters but
            # alpha
            def fit_test2(k, alp):
                return an.fit_smoothed_bPL(k, A=EM_mod[indkp], a=a, b=b,
                                           alp=alp, xc=ks[indkp])
            popt = opt.curve_fit(fit_test2, ks, EM_mod, p0=(alp))
            alp0=popt[0][0]
            if indt == 0:
                self.df['alp'] = alp0
                OmM_an = an.A_alpha(alp=alp0)*ks[indkp]*EM_mod[indkp]
                self.df['OmM (an)'] = OmM_an
                OmM_num = np.trapz(EM[indt, :], k)
                self.df['OmM (num)'] = OmM_num
                vA = np.sqrt(3/2*OmM_an)
                self.df['vA (an)'] = vA
                vA_num = np.sqrt(3/2*OmM_num)
                self.df['vA (num)'] = vA_num
                self.df['te (an)'] = 1/vA/ks[indkp]
                self.df['te (num)'] = 1/vA_num/ks[indkp]

    def magnetic_field_pars_t(self, a=4, b=5/3, kpeak0=1., Epeak0=1., alp=2.,
                              Nk=100000, fixb=True):

        import numpy as np
        import GW_analytical as an
        import scipy.optimize as opt

        #
        # use the smoothed broken power law models to fit the magnetic
        # spectra at each time step and obtain parameters by fixing all
        # parameters except amplitude and peak position (and in one of the
        # two cases, alpha)
        #

        # added option to allow for not fixing the slope b (sometimes it
        # is necessary to appropriate get peak position and amplitude)
        if fixb:
            def fit_test_fixalp(k, A, kpeak):
                return an.fit_smoothed_bPL(k, A=A, a=a, b=b, alp=alp, xc=kpeak)
            def fit_test(k, A, kpeak, alp):
                return an.fit_smoothed_bPL(k, A=A, a=a, b=b, alp=alp, xc=kpeak)
        else:
            def fit_test_fixalp(k, A, kpeak, b):
                return an.fit_smoothed_bPL(k, A=A, a=a, b=b, alp=alp, xc=kpeak)
            def fit_test(k, A, kpeak, alp, b):
                return an.fit_smoothed_bPL(k, A=A, a=a, b=b, alp=alp, xc=kpeak)

        k = self.spectra.get('k')
        tk = self.spectra.get('t_mag')
        if 'mag' in self.spectra_avail:
            EM = self.spectra.get('mag')
            EMpeaks = np.zeros(len(tk))
            kpeaks = np.zeros(len(tk))
            EMpeaks_fa = np.zeros(len(tk))
            kpeaks_fa = np.zeros(len(tk))
            if not fixb:
                bs = np.zeros(len(tk))
                bs_fa = np.zeros(len(tk))
            alphas = np.zeros(len(tk))
            ks = np.logspace(np.log10(k[1]), np.log10(k[-1]), Nk)
            EM_mod_fa = np.zeros((len(tk), len(ks)))
            EM_mod = np.zeros((len(tk), len(ks)))
            AA0_fa = Epeak0
            kp0_fa = kpeak0
            AA0 = AA0_fa
            kp0 = kp0_fa
            alp0 = alp
            b0_fa = b
            b0 = b
            # fit at each time step the amplitude and spectral peak for fixed
            # alpha (_fa) or allowing alpha to vary too
            for i in range(0, len(tk)):
                if fixb:
                    popt_fa = opt.curve_fit(fit_test_fixalp, k[1:], EM[i, 1:],
                                            p0=(AA0_fa, kp0_fa))
                else:
                    popt_fa = opt.curve_fit(fit_test_fixalp, k[1:], EM[i, 1:],
                                            p0=(AA0_fa, kp0_fa, b0_fa))
                EMpeaks_fa[i] = popt_fa[0][0]
                kpeaks_fa[i] = popt_fa[0][1]
                if not fixb: bs_fa[i] = popt_fa[0][2]
                AA0_fa = EMpeaks_fa[i]
                kp0_fa = kpeaks_fa[i]
                if not fixb: b0_fa = bs_fa[i]
                if fixb:
                    EM_mod_fa[i, :] = fit_test_fixalp(ks, EMpeaks_fa[i],
                                                      kpeaks_fa[i])
                else:
                    EM_mod_fa[i, :] = fit_test_fixalp(ks, EMpeaks_fa[i],
                                                      kpeaks_fa[i], bs_fa[i])

                if fixb:
                    popt = opt.curve_fit(fit_test, k[1:], EM[i, 1:],
                                         p0=(AA0, kp0, alp0))
                else:
                    popt = opt.curve_fit(fit_test, k[1:], EM[i, 1:],
                                         p0=(AA0, kp0, alp0, b0))

                EMpeaks[i] = popt[0][0]
                kpeaks[i] = popt[0][1]
                alphas[i] = popt[0][2]
                if not fixb: bs[i] = popt[0][3]
                AA0 = EMpeaks[i]
                kp0 = kpeaks[i]
                alp0 = alphas[i]
                if not fixb: b0 = bs[i]

                if fixb:
                    EM_mod[i, :] = fit_test(ks, EMpeaks[i], kpeaks[i],
                                            alphas[i])
                else:
                    EM_mod[i, :] = fit_test(ks, EMpeaks[i], kpeaks[i],
                                            alphas[i], bs[i])

            # compute decay exponents as function of time
            dte = self.df['te (num)'][0]
            q_exp = np.log(kpeaks[0]/kpeaks)/\
                            np.log(1 + (tk - self.tini)/dte)
            q_exp_fa = np.log(kpeaks_fa[0]/kpeaks_fa)/\
                            np.log(1 + (tk - self.tini)/dte)
            p_exp = np.log(EMpeaks[0]/EMpeaks)/\
                            np.log(1 + (tk - self.tini)/dte)
            p_exp_fa = np.log(EMpeaks_fa[0]/EMpeaks_fa)/\
                            np.log(1 + (tk - self.tini)/dte)

            if fixb:
                self.spectra.update({'EMpeaks_fa': EMpeaks_fa,
                                 'kpeaks_fa': kpeaks_fa,
                                 'EMpeaks': EMpeaks, 'kpeaks': kpeaks,
                                 'alphas': alphas, 'EM_model': EM_mod,
                                 'EM_model_fa': EM_mod_fa,
                                 'k_EM_model': ks, 'q_exp': q_exp,
                                 'q_exp_fa': q_exp_fa, 'p_exp': p_exp,
                                 'p_exp_fa': p_exp_fa})
            else:
                self.spectra.update({'bs': bs, 'bs_fa': bs_fa,
                                 'EMpeaks_fa_vb': EMpeaks_fa,
                                 'kpeaks_fa_vb': kpeaks_fa,
                                 'EMpeaks_vb': EMpeaks, 'kpeaks_vb': kpeaks,
                                 'alphas_vb': alphas, 'EM_model_vb': EM_mod,
                                 'EM_model_fa_vb': EM_mod_fa,
                                 'k_EM_model': ks, 'q_exp_vb': q_exp,
                                 'q_exp_fa_vb': q_exp_fa, 'p_exp_vb': p_exp,
                                 'p_exp_fa_vb': p_exp_fa})
            if not 'EM_model' in self.spectra_avail:
                self.spectra_avail.append('EM_model')
            if not 'EM_model_fa' in self.spectra_avail:
                self.spectra_avail.append('EM_model_fa')

    def compute_PiM_model(self, a=4, b=5/3, kpeak0=1., Epeak0=1., alp=2.,
                          Nk=100000, indf=-1, lowk=-4, highk=4):

        """
        Function that computes the spectrum of Pi at all times
        (only available for mag).

        Only available for alpha = 2 for now!!
        """

        import numpy as np
        import GW_analytical as an

        # Check if EM_model has been computed (and hence EMpeaks and
        # kpeaks), otherwise compute it
        if 'EM_model' not in self.spectra_avail:
            self.magnetic_field_pars_t(a=a, b=b, kpeak0=kpeak0,
                                       Epeak0=Epeak0, alp=alp,
                                       Nk=Nk)

        ks = self.spectra.get('k_EM_model')*1.
        tk = self.spectra.get('t_mag')
        # use values for a fixed alpha since PiM is not available for all alpha
        EMpeaks = self.spectra.get('EMpeaks_fa')
        kpeaks = self.spectra.get('kpeaks_fa')

        # check if we want to add k points out of the range for EM_model
        llow = False
        lhigh = False
        if lowk < np.log10(ks[0]):
            llow = True
            k_add_l = np.logspace(lowk, np.log10(ks[0]),
                                int(1000*(np.log10(ks[0]) - lowk)))
            ks = np.append(k_add_l, ks)
        if highk > np.log10(ks[-1]):
            lhigh=True
            k_add_h = np.logspace(np.log10(ks[-1]), highk,
                                int(1000*(highk - np.log10(ks[-1]))))
            ks = np.append(ks, k_add_h)

        # compute PiM from kpeak and EMpeak at each time
        PiMs = np.zeros((len(tk), len(ks)))
        if indf < 0: ltk = len(tk)
        else: ltk = int(indf)

        # compute PiM from normalized PiM (alpha = 2) and interpolate in ks
        for i in range(0, ltk):
            kf_Pi, PiM = an.compute_PiM(EMpeaks[i], kpeaks[i], alpha=2.,
                                        alp_str='2')
            PiMs[i, :] = np.interp(ks, kf_Pi, PiM)

        if not 'PiM_model_EM' in self.spectra_avail:
            self.spectra_avail.append('PiM_model_EM')
        self.spectra.update({'PiM_model_EM': PiMs, 'k_PiM_model_EM': ks})

def interpolate_ts(t1, t2, sp1, sp2):

    """
    Function that compares 2 series to compare them and interpolates both
    to the same time array. They should have same initial time.

    Arguments:
        t1, t2 -- time arrays of spectra 1 and 2
        sp1, sp2 -- spectrum arrays of spectra 1 and 2
    Returns:
        t_i -- common time array
        sp1, sp2 -- spectrum arrays of spectra 1 and 2 after interpolating
                    to have a common time array
    """

    import numpy as np

    # compare the end time of both runs
    end_t1 = t1[-1]
    end_t2 = t2[-1]
    if end_t1 < end_t2:
        sp2 = np.interp(t1, t2, sp2)
        t_i = t1*1
    else:
        sp1 = np.interp(t2, t1, sp1)
        t_i = t2*1

    return t_i, sp1, sp2
        
############################# NOT PUBLIC ON GITHUB ########################

def read_df_runs(dd='all_runs', dir=''):

    import os
    import pandas as pd

    if dir == '': dir = dir0

    CWD = os.getcwd()
    os.chdir(dir)
    df = pd.read_csv(dd + '.csv')
    os.chdir(CWD)

    return df

def update_df_runs(new_df, dd='all_runs', dir='', dd2='all_runs_new'):

    import os
    import pandas as pd

    if dir == '': dir = dir0

    df = read_df_runs(dd=dd, dir=dir)
    dfNew = df.merge(new_df, on=['project', 'dir', 'name'], left_index=True,
                     right_index=True,
                     how='outer', suffixes=('', '_y'))
    dfNew.drop(dfNew.filter(regex='_y$').columns, axis=1, inplace=True)
    CWD = os.getcwd()
    os.chdir(dir)
    dfNew.to_csv(dd2 + '.csv', index=False)
    os.chdir(CWD)

    return dfNew

def correct_ini1_3(runs):

    """
    Function that corrects the wrong normalization factor in runs ini1, ini2,
    and ini3.
    """

    print('\n')
    print('CORRECTING NORMALIZATION OF EGW IN RUNS ini1, ini2, AND ini3',
          'BY A FACTOR 32 pi/6 \n')

    rrs = ['ini1', 'ini2', 'ini3']
    for i in rrs:
        run = runs.get(i)
        EEGW = run.ts.get('EEGW')
        EEGW *= 16*np.pi/3
        run.ts.update({'EEGW':EEGW})
