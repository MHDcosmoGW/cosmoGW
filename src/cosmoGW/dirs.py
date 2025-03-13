"""
dirs.py is a Python routine that contains a function to create a dictionary
between the identifying names of runs in the literature with the corresponding
names of their directories, which should be included in the same directory.

The runs correspond to Pencil Code simulations of gravitational waves
in the early universe; see cosmoGW project (https://github.com/AlbertoRoper/cosmoGW)
for details.

Author: Alberto Roper Pol
Created: 01/01/2021
Updated: 01/11/2023 (new release of the cosmoGW code)

Main references are:

RPMBKK20 - A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili, A. Kosowsky,
"Numerical simulations of gravitational waves from early-universe turbulence,"
Phys. Rev. D 102, 083512 (2020), arXiv:1903.08585

KBGMRP21 - T. Kahniashvili, A. Brandenburg, G. Gogoberidze, S. Mandal, A. Roper Pol,
"Circular polarization of gravitational waves from early-universe helical
turbulence", Phys. Rev. Research 3, 013193 (2021), arxiv:2011.05556

RPMBK22 - A. Roper Pol, S. Mandal, A. Brandenburg, T. Kahniashvili,
"Polarization of gravitational waves from helical MHD turbulent sources,"
JCAP 04 (2022), 019, arXiv:2107.05356

HRPB21 - Y. He, A. Roper Pol, A. Brandenburg, "Leading-order nonlinear
gravitational waves from reheating magnetogeneses", arxiv:2110.14456 (2021)

RPCNS22 - A. Roper Pol, C. Caprini, A. Neronov, D. Semikoz, "The gravitational wave
signal from primordial magnetic fields in the Pulsar Timing Array frequency band,"
Phys. Rev. D 105, 123502 (2022), arXiv:2201.05630

HRPB23 - Y. He, A. Roper Pol, A. Brandenburg, "Modified propagation of gravitational
waves from the early radiation era," JCAP 06 (2023), 025, arXiv:2212.06082

"""

def read_dirs(proj, dirs={}, ext=''):

    """
    Function that returns an updated dictionary of the run directories.

    Arguments:
        proj -- project at which the directories correspond
        dirs -- initial directory with the runs of the projects to be updated
                (default is empty directory)
        ext -- extension selecting specific subsets of a project (by default
               reads all runs of a project)

    The available projects proj are:

        'PRD_1903_08585'  -- runs of RPMBKK20

                              The subsets of runs are: 'ini' (initially given),
                              'hel' (helical forced), 'noh' (non-helical forced),
                              and 'ac' (acoustic turbulence)

        'PRR_2011_05556'  -- runs of KBGMRP21

                             The subsets of runs are: 'K' (kinetic), 'M' (magnetic), and
                             'nohel'

        'JCAP_2107_05356' -- runs of RPMBK22

                             The subsets of runs are: 'ini' (initially given) and
                             'forc' (initially driven).

        'PRD_2110_14456'  -- runs of HRPB21

                             The subsets of runs are: 'helical_b27',
                             'helical_b73', 'nonhelical_b27',
                             'nonhelical_b73', and 'helical_b17'
                             for helical and non-helical runs with
                             beta = 1.7, 2.7, and 7.3, and '_toff'
                             for extended runs (up to t = 10).

        'PRD_2201_05630'  -- runs of RPCNS22

        'JCAP_2212_06082' -- runs of HRPB23

                             The subsets of runs are: 'M0', 'M1', 'M2', 'M3' corresponding
                             to each of the 0-4 parameterizations of alpM with time.
                             To choose runs with large domain size, specify '_lowk' after the
                             corresponding subset of runs.

    Returns:
        dirs -- updated dictionary of directories
    """

    ########################## PRD_1903_08585 ##################################
    
    '''
    Reference for PRD_1903_08585 is RPMBKK20
    
    Zenodo: https://zenodo.org/record/3692072
    '''

    if proj=='PRD_1903_08585':

        if ext=='ini' or ext=='':

            dirs.update({'ini1': 'M1152e_exp6k4_M4b'})
            dirs.update({'ini2': 'M1152e_exp6k4'})
            dirs.update({'ini3': 'M1152e_exp6k4_k60b'})

        if ext=='hel' or ext=='':

            dirs.update({'hel1': 'F1152d2_sig1_t11_M2c_double'})
            dirs.update({'hel2': 'F1152a_sig1_t11d_double'})
            dirs.update({'hel3': 'F1152a_sig1'})
            dirs.update({'hel4': 'F1152a_k10_sig1'})

        if ext=='noh' or ext=='':

            dirs.update({'noh1': 'F1152b_sig0_t11_M4'})
            dirs.update({'noh2': 'F1152a_sig0_t11b'})

        if ext=='ac' or ext=='':

            dirs.update({'ac1':'E1152e_t11_M4d_double'})
            dirs.update({'ac2':'E1152e_t11_M4a_double'})
            dirs.update({'ac3':'E1152e_t11_M4e_double'})

    ########################## PRR_2011_05556 ##################################
    
    '''
    Reference for PRR_2011_05556 is KBGMRP21
    
    Zenodo: https://zenodo.org/record/4256906
    '''

    if proj=='PRR_2011_05556':

        if ext=='K' or ext=='':

            dirs.update({'K0': 'K512sig0_k6_ramp1a'})
            dirs.update({'K01_c': 'K512sig01_k6_ramp1c'})
            dirs.update({'K01_a': 'K512sig01_k6_ramp1a'})
            dirs.update({'K03': 'K512sig03_k6_ramp1a'})
            dirs.update({'K05': 'K512sig05_k6_ramp1a'})
            dirs.update({'K1': 'K512sig1_k6_ramp1a'})

        if ext=='M' or ext=='':

            dirs.update({'M0': 'M512sig0_k6_ramp1a'})
            dirs.update({'M01_c': 'M512sig01_k6_ramp1c'})
            dirs.update({'M01_b': 'M512sig01_k6_ramp1b'})
            dirs.update({'M03': 'M512sig03_k6_ramp1a'})
            dirs.update({'M05': 'M512sig05_k6_ramp1a'})
            dirs.update({'M1': 'M512sig1_k6_ramp1a'})

        if ext=='nohel' or ext=='':

            dirs.update({'nohel_tau01': 'F1152a_sig0_t11_M4_ramp01b'})
            dirs.update({'nohel_tau02': 'F1152a_sig0_t11_M4_ramp02a'})
            dirs.update({'nohel_tau05': 'F1152a_sig0_t11_M4_ramp05a'})
            dirs.update({'nohel_tau1': 'F1152a_sig0_t11_M4_ramp1a'})
            dirs.update({'nohel_tau2': 'F1152a_sig0_t11_M4_ramp2a'})

    ########################## JCAP_2107_05356 ################################
    
    '''
    Reference for JCAP_2107_05356 is RPMBK22
    
    Zenodo: https://zenodo.org/record/5525504
    '''

    if proj=='JCAP_2107_05356':

        if ext=='ini' or ext=='':

            dirs.update({'i_s01': 'M1152e_exp6k4_sig01'})
            dirs.update({'i_s03': 'M1152e_exp6k4_sig03'})
            dirs.update({'i_s05': 'M1152e_exp6k4_sig05'})
            dirs.update({'i_s07': 'M1152e_exp6k4_sig07'})
            dirs.update({'i_s1': 'M1152e_exp6k4'})

        if ext=='forc' or ext=='':

            dirs.update({'f_s001_neg': 'F1152sigm001a'})
            dirs.update({'f_s001': 'F1152sig001a'})
            dirs.update({'f_s03': 'F1152sig03a'})
            dirs.update({'f_s05': 'F1152sig05a'})
            dirs.update({'f_s07': 'F1152sig07c'})
            dirs.update({'f_s1_neg': 'F1152sigm1a'})

    ########################### PRD_2110_14456 ################################
    
    '''
    Reference for PRD_2110_14456 is HRPB21
    
    Zenodo: https://zenodo.org/record/5603013
    '''

    if proj=='PRD_2110_14456':

        if 'nonhelical_b73' in ext or ext=='':

            if 'toff' in ext or ext=='':

                dirs.update({'A1_nl_toff': 'P512b73c_nlin2_nhel_e002new_toff'})
                dirs.update({'A1_l_toff': 'P512b73c_lin2_nhel_e002new_toff'})
                dirs.update({'A2_nl_toff': 'P512b73c_nlin2_nhel_e01new_toff'})
                dirs.update({'A2_l_toff': 'P512b73c_lin2_nhel_e01new_toff'})

            elif 'toff' not in ext or ext=='':

                dirs.update({'A1_nl': 'P512b73c_nlin2_nhel_e002new'})
                dirs.update({'A1_l': 'P512b73c_lin2_nhel_e002new'})
                dirs.update({'A2_nl': 'P512b73c_nlin2_nhel_e01new'})
                dirs.update({'A2_l': 'P512b73c_lin2_nhel_e01new'})
                dirs.update({'A3_nl': 'P512b73c_nlin2_nhel_e1new'})
                dirs.update({'A3_l': 'P512b73c_lin2_nhel_e1new'})
                dirs.update({'A4_nl': 'P512b73c_nlin2_nhel_e10new'})
                dirs.update({'A4_l': 'P512b73c_lin2_nhel_e10new'})

        if 'nonhelical_b27' in ext or ext=='':

            if 'toff' in ext or ext=='':

                dirs.update({'B1_nl_toff': 'P512b27c_nlin2_nhel_e002new_toff'})
                dirs.update({'B1_l_toff': 'P512b27c_lin2_nhel_e002new_toff'})
                dirs.update({'B2_nl_toff': 'P512b27c_nlin2_nhel_e01new_toff'})
                dirs.update({'B2_l_toff': 'P512b27c_lin2_nhel_e01new_toff'})

            elif 'toff' not in ext or ext=='':

                dirs.update({'B1_nl': 'P512b27c_nlin2_nhel_e002new'})
                dirs.update({'B1_l': 'P512b27c_lin2_nhel_e002new'})
                dirs.update({'B2_nl': 'P512b27c_nlin2_nhel_e01new'})
                dirs.update({'B2_l': 'P512b27c_lin2_nhel_e01new'})
                dirs.update({'B3_nl': 'P512b27c_nlin2_nhel_e1new'})
                dirs.update({'B3_l': 'P512b27c_lin2_nhel_e1new'})
                dirs.update({'B4_nl': 'P512b27c_nlin2_nhel_e10new'})
                dirs.update({'B4_l': 'P512b27c_lin2_nhel_e10new'})

        if 'helical_b73' in ext or ext=='':

            if 'toff' in ext or ext=='':

                dirs.update({'C1_nl_toff': 'P512b73c_nlin2_e002new_toff'})
                dirs.update({'C1_l_toff': 'P512b73c_lin2_e002new_toff'})
                dirs.update({'C2_nl_toff': 'P512b73c_nlin2_e01new_toff'})
                dirs.update({'C2_l_toff': 'P512b73c_lin2_e01new_toff'})

            elif 'toff' not in ext or ext=='':

                dirs.update({'C1_nl': 'P512b73c_nlin2_e002new'})
                dirs.update({'C1_l': 'P512b73c_lin2_e002new'})
                dirs.update({'C2_nl': 'P512b73c_nlin2_e01new'})
                dirs.update({'C2_l': 'P512b73c_lin2_e01new'})
                dirs.update({'C3_nl': 'P512b73c_nlin2_e1new'})
                dirs.update({'C3_l': 'P512b73c_lin2_e1new'})
                dirs.update({'C4_nl': 'P512b73c_nlin2_e10new'})
                dirs.update({'C4_l': 'P512b73c_lin2_e10new'})

        if 'helical_b27' in ext or ext=='':

            if 'toff' in ext or ext=='':

                dirs.update({'D1_nl_toff': 'P512b27c_nlin2_e002new_toff'})
                dirs.update({'D1_l_toff': 'P512b27c_lin2_e002new_toff'})
                dirs.update({'D2_nl_toff': 'P512b27c_nlin2_e01new_toff'})
                dirs.update({'D2_l_toff': 'P512b27c_lin2_e01new_toff'})

            elif 'toff' not in ext or ext=='':

                dirs.update({'D1_nl': 'P512b27c_nlin2_e002new'})
                dirs.update({'D1_l': 'P512b27c_lin2_e002new'})
                dirs.update({'D2_nl': 'P512b27c_nlin2_e01new'})
                dirs.update({'D2_l': 'P512b27c_lin2_e01new'})
                dirs.update({'D3_nl': 'P512b27c_nlin2_e1new'})
                dirs.update({'D3_l': 'P512b27c_lin2_e1new'})
                dirs.update({'D4_nl': 'P512b27c_nlin2_e10new'})
                dirs.update({'D4_l': 'P512b27c_lin2_e10new'})

        if 'helical_b17' in ext or ext=='':

            if 'toff' in ext or ext=='':

                dirs.update({'E1_nl_toff': 'P512b17c_k02b_nlin2_e002new_toff'})
                dirs.update({'E1_l_toff': 'P512b17c_k02b_lin2_e002new_toff'})
                dirs.update({'E2_nl_toff': 'P512b17c_k02b_nlin2_e01new_toff'})
                dirs.update({'E2_l_toff': 'P512b17c_k02b_lin2_e01new_toff'})

            elif 'toff' not in ext or ext=='':

                dirs.update({'E1_nl': 'P512b17c_k02b_nlin2_e002new'})
                dirs.update({'E1_l': 'P512b17c_k02b_lin2_e002new'})
                dirs.update({'E2_nl': 'P512b17c_k02b_nlin2_e01new'})
                dirs.update({'E2_l': 'P512b17c_k02b_lin2_e01new'})
                dirs.update({'E3_nl': 'P512b17c_k02b_nlin2_e1new'})
                dirs.update({'E3_l': 'P512b17c_k02b_lin2_e1new'})

    ########################### PRD_2201_05630 ################################
    
    '''
    Reference for PRD_2201_05630 is RPCNS22
    
    Zenodo: https://zenodo.org/record/5782752
    '''

    if proj == 'PRD_2201_05630':

        dirs.update({'A1': 'M768_k15_sig0_OmM01_k0_033'})
        dirs.update({'A2': 'M768_k15_sig0_OmM01_k0_017'})
        dirs.update({'B': 'M768_k11_sig0_OmM01_k0_033'})
        dirs.update({'C1': 'M768_k8_sig0_OmM01_k0_033'})
        dirs.update({'C2': 'M768_k8_sig0_OmM01_k0_017'})
        dirs.update({'D1': 'M768_k7_sig0_OmM01_k0_033'})
        dirs.update({'D2': 'M768_k7_sig0_OmM01_k0_017'})
        dirs.update({'E1': 'M512_k6_sig0_OmM8em3_k0_05'})
        dirs.update({'E2': 'M512_k6_sig0_OmM8em3_k0_02'})
        dirs.update({'E3': 'M512_k6_sig0_OmM8em3_k0_01'})
        dirs.update({'E4': 'M512_k6_sig0_OmM8em3_k0_0075'})
        dirs.update({'E5': 'M512_k6_sig0_OmM8em3_k0_0033'})
        
    ########################### JCAP_2212_06082 ################################
    
    '''
    Reference for JCAP_2212_06082 is HRPB23
    
    Zenodo: https://zenodo.org/record/7408601
    '''

    if proj == 'JCAP_2212_06082':

        if 'M0' in ext or ext=='':
            
            if 'lowk' not in ext or ext=='':

                dirs.update({'M0A' : '46000cosmo_1D_alpM_m05'})
                dirs.update({'M0A_LD2' : '46000cosmo_1D_alpM_m05_LD2'})
                dirs.update({'M0B' : '46000cosmo_1D_alpM_m03'})
                dirs.update({'M0B_LD2' : '46000cosmo_1D_alpM_m03_LD2'})
                dirs.update({'M0C' : '46000cosmo_1D_alpM_m01'})
                dirs.update({'M0C_LD2' : '46000cosmo_1D_alpM_m01_LD2'})
                dirs.update({'M0D' : '46000cosmo_1D_alpM_m001'})
                dirs.update({'M0D_LD2' : '46000cosmo_1D_alpM_m001_LD2'})
                dirs.update({'M0E' : '46000cosmo_1D_alpM_p01'})
                dirs.update({'M0E_LD2' : '46000cosmo_1D_alpM_p01_LD2'})
                dirs.update({'M0F' : '46000cosmo_1D_alpM_p03'})
                dirs.update({'M0F_LD2' : '46000cosmo_1D_alpM_p03_LD2'})

            if 'lowk' in ext or ext=='':
                
                dirs.update({'M0A_lowk' : '46000cosmo_1D_alpM_m05_wav1e7'})
                dirs.update({'M0A_lowk_LD2' : '46000cosmo_1D_alpM_m05_wav1e7_LD2'})
                dirs.update({'M0B_lowk' : '46000cosmo_1D_alpM_m03_wav1e7'})
                dirs.update({'M0B_lowk_LD2' : '46000cosmo_1D_alpM_m03_wav1e7_LD2'})
                dirs.update({'M0C_lowk' : '46000cosmo_1D_alpM_m01_wav1e7'})
                dirs.update({'M0C_lowk_LD2' : '46000cosmo_1D_alpM_m01_wav1e7_LD2'})
                dirs.update({'M0D_lowk' : '46000cosmo_1D_alpM_m001_wav1e7'})
                dirs.update({'M0D_lowk_LD2' : '46000cosmo_1D_alpM_m001_wav1e7_LD2'})
                dirs.update({'M0E_lowk' : '46000cosmo_1D_alpM_p01_wav1e7'})
                dirs.update({'M0E_lowk_LD2' : '46000cosmo_1D_alpM_p01_wav1e7_LD2'})
                dirs.update({'M0F_lowk' : '46000cosmo_1D_alpM_p03_wav1e7'})
                dirs.update({'M0F_lowk_LD2' : '46000cosmo_1D_alpM_p03_wav1e7_LD2'})

        if 'M1' in ext or ext=='':

            dirs.update({'M1A' : '46000cosmo_1D_alpM_m05_a2'})
            dirs.update({'M1A_LD2' : '46000cosmo_1D_alpM_m05_a2_LD2'})
            dirs.update({'M1B' : '46000cosmo_1D_alpM_m03_a2'})
            dirs.update({'M1B_LD2' : '46000cosmo_1D_alpM_m03_a2_LD2'})
            dirs.update({'M1C' : '46000cosmo_1D_alpM_m01_a2'})
            dirs.update({'M1C_LD2' : '46000cosmo_1D_alpM_m01_a2_LD2'})
            dirs.update({'M1D' : '46000cosmo_1D_alpM_p01_a04'})
            dirs.update({'M1D_LD2' : '46000cosmo_1D_alpM_p01_a04_LD2'})
            dirs.update({'M1E' : '46000cosmo_1D_alpM_p03_a04'})
            dirs.update({'M1E_LD2' : '46000cosmo_1D_alpM_p03_a04_LD2'})

        if 'M2' in ext or ext=='':

            dirs.update({'M2A' : '46000cosmo_1D_alpM_m05_Lambda'})
            dirs.update({'M2A_LD2' : '46000cosmo_1D_alpM_m05_Lambda_LD2'})
            dirs.update({'M2B' : '46000cosmo_1D_alpM_m03_Lambda'})
            dirs.update({'M2B_LD2' : '46000cosmo_1D_alpM_m03_Lambda_LD2'})
            dirs.update({'M2C' : '46000cosmo_1D_alpM_m01_Lambda'})
            dirs.update({'M2C_LD2' : '46000cosmo_1D_alpM_m01_Lambda_LD2'})
            dirs.update({'M2D' : '46000cosmo_1D_alpM_p01_Lambda'})
            dirs.update({'M2D_LD2' : '46000cosmo_1D_alpM_p01_Lambda_LD2'})
            dirs.update({'M2E' : '46000cosmo_1D_alpM_p03_Lambda'})
            dirs.update({'M2E_LD2' : '46000cosmo_1D_alpM_p03_Lambda_LD2'})

        if 'M3' in ext or ext=='':
            
            if 'lowk' not in ext or ext=='':

                dirs.update({'M3A' : '46000cosmo_1D_alpM_m05_mat'})
                dirs.update({'M3A_LD2' : '46000cosmo_1D_alpM_m05_mat_LD2'})
                dirs.update({'M3B' : '46000cosmo_1D_alpM_m03_mat'})
                dirs.update({'M3B_LD2' : '46000cosmo_1D_alpM_m03_mat_LD2'})
                dirs.update({'M3C' : '46000cosmo_1D_alpM_m01_mat'})
                dirs.update({'M3C_LD2' : '46000cosmo_1D_alpM_m01_mat_LD2'})
                dirs.update({'M3D' : '46000cosmo_1D_alpM_m001_mat'})
                dirs.update({'M3D_LD2' : '46000cosmo_1D_alpM_m001_mat_LD2'})
                dirs.update({'M3E' : '46000cosmo_1D_alpM_p01_mat'})
                dirs.update({'M3E_LD2' : '46000cosmo_1D_alpM_p01_mat_LD2'})
                dirs.update({'M3F' : '46000cosmo_1D_alpM_p03_mat'})
                dirs.update({'M3F_LD2' : '46000cosmo_1D_alpM_p03_mat_LD2'})

            if 'lowk' in ext or ext=='':
                
                dirs.update({'M3A_lowk' : '46000cosmo_1D_alpM_m05_mat_wav1e7'})
                dirs.update({'M3A_lowk_LD2' : '46000cosmo_1D_alpM_m05_mat_wav1e7_LD2'})
                dirs.update({'M3B_lowk' : '46000cosmo_1D_alpM_m03_mat_wav1e7'})
                dirs.update({'M3B_lowk_LD2' : '46000cosmo_1D_alpM_m03_mat_wav1e7_LD2'})
                dirs.update({'M3C_lowk' : '46000cosmo_1D_alpM_m01_mat_wav1e7'})
                dirs.update({'M3C_lowk_LD2' : '46000cosmo_1D_alpM_m01_mat_wav1e7_LD2'})
                dirs.update({'M3D_lowk' : '46000cosmo_1D_alpM_m001_mat_wav1e7'})
                dirs.update({'M3D_lowk_LD2' : '46000cosmo_1D_alpM_m001_mat_wav1e7_LD2'})
                dirs.update({'M3E_lowk' : '46000cosmo_1D_alpM_p01_mat_wav1e7'})
                dirs.update({'M3E_lowk_LD2' : '46000cosmo_1D_alpM_p01_mat_wav1e7_LD2'})
                dirs.update({'M3F_lowk' : '46000cosmo_1D_alpM_p03_mat_wav1e7'})
                dirs.update({'M3F_lowk_LD2' : '46000cosmo_1D_alpM_p03_mat_wav1e7_LD2'})

    ########################## unpublished directories ########################

    if proj == 'horndeski':
        dirs.update({'A': '64000_cosmo_GR_k20_kw2Em3_40proc'})
        dirs.update({'A2': '64000_cosmo_GR_k20_kw2Em3_40proc_smallk'})
        dirs.update({'B': '64000_cosmo_GR_k20_kw2Em3'})
        dirs.update({'D': '64000_cosmo_modGR_alpMm001_k20_kw2Em3_40proc'})
        dirs.update({'D2': '64000_cosmo_modGR_alpMm001_k20_kw2Em3_40proc_smallk'})
        dirs.update({'D3': '64000_cosmo_modGR_alpMm001_k20_kw2Em3_40proc_smallk2'})
        dirs.update({'D4': '64000_cosmo_modGR_alpMm001_k20_kw2Em3_40proc_smallk3'})
        dirs.update({'D5': '64000_cosmo_modGR_alpMm001_k20_kw2Em3_40proc_shortt'})
        dirs.update({'E': '64000_cosmo_modGR_alpMm005_k20_kw2Em3_40proc'})
        dirs.update({'C': '64000_cosmo_modGR_alpMm01_k20_kw2Em3_40proc'})
        dirs.update({'F': '64000_cosmo_modGR_alpTm05_k20_kw2Em3_40proc'})

    if proj == 'nanograv':

        dirs.update({'run1': 'M256k10_sig0_OmM0003'})
        dirs.update({'run2': 'M256k10_sig0_OmM09'})
        dirs.update({'run3': 'M256k60_sig0_OmM00004'})
        dirs.update({'run4': 'M256k60_sig0_OmM00004_smallt'})
        dirs.update({'run5': 'M256k10_sig0_OmM03'})
        dirs.update({'run6': 'M256k10_sig0_OmM003'})
        dirs.update({'run7':'M256k10_sig0_OmM003_smallt'})
        dirs.update({'run8': 'M256k10_sig0_OmM02'})
        dirs.update({'run9': 'M256k10_sig0_OmM02_large_visc'})
        dirs.update({'run11': 'M256k10_sig0_OmM01'})
        dirs.update({'run12': 'M256k10_sig0_OmM01b'})
        dirs.update({'run13_sd': 'M512k10_sig0_OmM0007_smalldom'})
        dirs.update({'run13_ld': 'M512k10_sig0_OmM0007_largedom'})
        dirs.update({'run13': 'M512k10_sig0_OmM0007_upd'})
        dirs.update({'run14': 'M512_k6_sig0_OmM001'})
        dirs.update({'run14_08': 'M512_k6_sig0_OmM0008'})
        dirs.update({'run14_k05': 'M512_k6_sig0_OmM0008b'})
        dirs.update({'run14_k02': 'M512_k6_sig0_OmM0008c'})
        dirs.update({'run14_k01': 'M512_k6_sig0_OmM0008d'})
        dirs.update({'run14_k0075': 'M512_k6_sig0_OmM0008e'})
        dirs.update({'run14_k0033': 'M512_k6_sig0_OmM0008f'})
        dirs.update({'run15': 'M512_k30_sig0_OmM00045'})
        dirs.update({'run16': 'M768_k60_sig0_OmM00035'})
        dirs.update({'run17': 'M768_k90_sig0_OmM00035'})
        dirs.update({'run18': 'M768_k14_sig0_OmM01'})
        dirs.update({'run18b': 'M768_k14_sig0_OmM01_ft8_k033'})
        dirs.update({'run18c': 'M768_k14_sig0_OmM01_k017'})
        dirs.update({'run19': 'M768_k6_sig0_OmM01'})
        dirs.update({'run19b': 'M768_k6_sig0_OmM01_k016'})
        dirs.update({'run20': 'M768_k10_sig0_OmM01'})
        dirs.update({'run21': 'M768_k8_sig0_OmM01'})
        dirs.update({'run21b': 'M768_k8_sig0_OmM01_k016'})

    if proj == 'GW_QCD_Clarke_2021_all':

        dirs.update({'QCD_run1': 'K512irro_k6_wav03_7p0em1d_512procs'})
        dirs.update({'QCD_run2': 'M512sig0_k6_wav3_2p3em2'})
        dirs.update({'QCD_run3': 'M1024sig0_k6_wav03_1p9em1'})
        dirs.update({'QCD_run4': 'M512sig0_k6_wav03_1p9em1'})
        dirs.update({'QCD_run5': 'M512sig1_k6_wav3_2p3em2'})
        dirs.update({'QCD_run6': 'M512sig1_k6_wav03_3p0em1'})
        dirs.update({'QCD_run7': 'M1024sig0_k6_wav1_1em1c_1024proc'})
        dirs.update({'QCD_run8': 'M512sig0_k6_wav03_3p0em1'})
        dirs.update({'QCD_run9': 'M512sig0_k6_wav1_2p0em2b'})
        dirs.update({'QCD_run10': 'M512sig1_k6_wav03_1p9em1b'})
        dirs.update({'QCD_run11': 'M512sig1_k6_wav1_1p0em1'})
        dirs.update({'QCD_run12': 'M512sig0_k6_wav03_1p0em1'})
        dirs.update({'QCD_run13': 'M512sig0_k6_wav10_1p0em2'})
        dirs.update({'QCD_run14': 'K512sig1_k6_wav03_4p0em1f'})
        dirs.update({'QCD_run15': 'M512sig1_k6_wav1b'})
        dirs.update({'QCD_run16': 'M512sig0_k6_wav1b'})
        dirs.update({'QCD_run17': 'M512sig1_k6_wav03_1p0em1'})
        dirs.update({'QCD_run18': 'M512sig1_k6_wav1_2p0em2'})
        dirs.update({'QCD_run19': 'M512sig1_k6_wav10_1p0em2b'})
        
    if proj == 'higgsless':
        
        dirs.update({'QCD_run1': 'K512irro_k6_wav03_7p0em1d_512procs'})

    return dirs
