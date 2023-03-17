"""
==============================================================
            gcm_toolkit Cloud functionalities
==============================================================
 This file contains multiple functionalities to deal with
 clouds in GCMs. Most of them are universal functionalites,
 some are GCM specific (see comments of individual functions).
==============================================================
"""
import os
import numpy as np
from glob import glob
from scipy.optimize import minimize
from ..core import writer as wrt


def cloud_opacities(wavelengths, cloud_radius, cloud_abundances, cloud_particle_density, volume_fraction):
    """
    Calculate cloud particle opacities using Mie-Theory. This function assumes spherical
    particles. The output are the absorption efficiency, scattering efficiency, and
    cross-section. One can calculate the Kappa values by multiplying the efficiencies
    with the cross-section.

    Parameters
    ----------
    wavelengths: np.ndarray[w]
        wavelengths to evaluate the opacities at [m].
    cloud_radius: np.ndarray[h]
        Mean cloud particle radius [m].
    cloud_abundances: np.ndarray[h]
        Mass fraction of clouds compared to the gas.
    cloud_particle_density: np.ndarray[h]
        Bulk density of cloud particles [kg/m^3]
        (for homogenous cloud particles, this is the same as the density of
        the cloud particle material)
    volume_fraction: dict[Species: np.ndarray[h]]
        Volume fractions of cloud species. All species included must be
        supported by the cloud_nk_data function.

    Returns
    -------
    Qabs : np.ndarray[h, w]
        Absorption coefficients
    Qsca : np.ndarray[h, w]
        Scattering coefficients
    c_sec : np.ndarray[h]
        cross sections
    """

    # check input
    if not(len(cloud_radius) == len(cloud_abundances) == len(cloud_particle_density)):
        raise ValueError('All cloud particle properties need to be '
                         'defined on the same pressure gird and have '
                         'the same number of entries. Given:'
                         '\n  len(cloud_radius) = ' + str(len(cloud_radius)) +
                         '\n  len(cloud_abundances) = ' + str(len(cloud_abundances)) +
                         '\n  len(cloud_particle_density) = ' + str(len(cloud_particle_density)))

    # low number handling
    cloud_particle_density[cloud_particle_density < 0] = 0
    cloud_radius[cloud_radius < 0] = 0

    # get array lengths
    len_h = len(cloud_radius)
    len_w = len(wavelengths)

    # prepare output arrays
    qabs = np.zeros((len_h, len_w))
    qsca = np.zeros((len_h, len_w))

    # calculate effective medium theory for each height and each wavelength
    eff = eff_func(volume_fraction, wavelengths)

    # load mie calculation library
    import miepython

    # loop over all inputs:
    for w in range(len_w):
        for h in range(len_h):
            print('mie', w, h)
            # prepare input parameters for mie calculation
            m = complex(eff[w, h, 0], -eff[w, h, 1])
            x = 2*np.pi*cloud_radius[h]/wavelengths[w]

            # calculate mie theory
            qe_temp, qs_temp, qback, g = miepython.mie(m, x)

            # save data
            qabs[h, w] = qe_temp - qs_temp
            qsca[h, w] = qs_temp

    # calculate cloud particle cross section
    c_sec = cloud_abundances * 3/4 / (cloud_radius * cloud_particle_density)

    # all done
    return qabs, qsca, c_sec


def eff_func(volume_fraction, wavelength):
    """
    Get the effective real and imaginary refractory index of a mixed
    cloud particle grain.

    Parameters
    ----------
    volume_fraction: dict[Species: np.ndarray[h]]
        Volume fractions of cloud species. All species included must be
        supported by the cloud_nk_data function.
    wavelength: np.ndarray[w]
        The wavelength at which the refractory index should be given [m].

    Returns
    -------
    ref_index: np.ndarray[w, h, 2]
        The refractory index (real, imaginary) with the first index
        being the wavelength, second index being the height point.
    """

    # get number of height points
    len_h = len(volume_fraction[list(volume_fraction.keys())[0]])
    # get number of wavelengths
    len_w = len(wavelength)
    # get number of bulk material species
    len_b = len(volume_fraction.keys())

    # prepare output
    eff = np.zeros((len_w, len_h, 2))

    # get all ref indexes
    work = np.zeros((len_b, len_h, len_w, 3))
    for k, key in enumerate(volume_fraction):
        temp = _cloud_nk_data(key, wavelength)
        for h in range(len_h):
            work[k, h, :, 1:] = temp

    # get all volume fractions
    for w in range(len_w):
        for h in range(len_h):
            for k, key in enumerate(volume_fraction):
                work[k, h, w, 0] = volume_fraction[key][h]

    # loop over all height points
    for w in range(len_w):
        for h in range(len_h):
            print('eff', w, h)
            # initial guess using linear approximation
            eff_0 = np.zeros((2,))
            eff_0[0] = sum(work[:, h, w, 0] * work[:, h, w, 1])
            eff_0[1] = sum(work[:, h, w, 0] * work[:, h, w, 2])

            # calculate effective medium theory with Bruggemann minimization
            eff[w, h] = minimize(_func, eff_0, args=work[:, h, w, :]).x

            # # LLL method
            # val = complex(0,0)
            # for k, key in enumerate(volume_fraction):
            #     val += work[k, 0] * complex(work[k, 1], work[k, 2])**(1/3)
            # val = val**3
            # tmp = np.sqrt(val.real**2 + val.imag**2)
            # eff[w, h, 0] = np.sqrt(0.5*(val.real + tmp))
            # eff[w, h, 1] = np.sqrt(0.5*(-val.real + tmp)) * val.imag/np.abs(val.imag)

    return eff


# define Bruggeman mixing rule
def _func(eff, work):
    """
    Bruggeman mixing function.

    Parameters
    ----------
    eff: np.ndarray[2]
        effective reffrective index. First index real part, second index
        imaginary part
    work: np.ndarray[s, 3]
        The volume mixing fraction, real refractive index and imaginary
        refractive index for each bulk material species s.

    Returns
    -------
    quality_of_fit : float
        A quality of fit parameter. Should be as close to zero as possible.
    """
    val = np.zeros((2,))
    for i, _ in enumerate(work):
        # transform from n,k to e describtion (a,c = real part, b,d = imaginary part)
        a = work[i, 1]**2 - work[i, 2]**2
        b = 2*work[i, 1]*work[i, 2]
        c = eff[0]
        d = eff[1]

        # calcualte vol_frac*(e - eff)/(e + 2*eff) for complex numbers
        temp = a**2 + 4*a*c + b**2 + 4*b*d + 4*(c**2 + d**2)
        val[0] += work[i, 0] * (a**2 + a*c + b**2 + b*d - 2*(c**2 + d**2)) / temp
        val[1] += work[i, 0] * (3*b*c - 3*a*d)/temp

    return val[0]**2 + val[1]**2


def _cloud_nk_data(species, wavelength):
    """
    Get the real (n) and imaginary (k) refractory index of a given species.

    Parameters
    ----------
    species: str
        Name of the cloud particle species. Currently supported:
        TiO2[s], Mg2SiO4[s], SiO[s], SiO2[s], Fe[s], Al2O3[s], CaTiO3[s], FeO[s],
        FeS[s], Fe2O3[s], MgO[s], MgSiO3[s], CaSiO3[s], Fe2SiO4[s], C[s], KCl[s]
    wavelength: Union[np.ndarray, float]
        The wavelength at which the refractory index should be given [m].

    Returns
    -------
    ref_index: np.ndarray[n, 2]
        The refractory index with the first index being the wavelenght and the second being n, k.
    """

    # gather all data files
    files = glob(os.path.dirname(__file__) + '/../nk_data/*.dat')

    # convert wavelengths to micron
    wavelength *= 1e6

    # prepare output
    if type(wavelength) is not np.ndarray:
        wavelength = np.asarray([wavelength])
    ref_index = np.zeros((len(wavelength), 2))

    # find species in files
    for fil in files:
        if species in fil:

            # get data
            content = open(fil, 'r').readlines()
            header = content[0].split()[1]
            data = content[5:]

            # loop over all wavelengths:
            for w, wave in enumerate(wavelength):
                # if desired wavelength is smaller than data, use the smallest
                # wavelength data available.
                if wave < float(data[0].split()[0]):
                    ref_index[w, 0] = float(data[0].split()[1])
                    ref_index[w, 1] = float(data[0].split()[2])
                    continue

                # if wavelength is within range, log-log interpolation
                for d, _ in enumerate(data):
                    cur_wave = float(data[d].split()[0])  # current wavelength
                    if wave < cur_wave:
                        nl = float(data[d-1].split()[1])  # lower n value
                        nh = float(data[d].split()[1])  # higher n value
                        kl = float(data[d-1].split()[2])  # lower k value
                        kh = float(data[d].split()[2])  # higher k value
                        prev_wave = float(data[d-1].split()[0])  # previous wavelength

                        # calculate interpolation
                        fac = np.log(wave/prev_wave)/np.log(cur_wave/prev_wave)
                        ref_index[w, 0] = np.exp(np.log(nl) + fac * np.log(nh/nl))
                        if kl <= 0 or kh < 0:
                            ref_index[w, 1] = 0
                        else:
                            ref_index[w, 1] = np.exp(np.log(kl) + fac * np.log(kh/kl))
                        break
                else:
                    # if wavelength is outside of range, extrapolate
                    if header == '.True.':
                        # conducting interpolation, log-log interpolation
                        # find 70% wavelength
                        cur_wave = float(data[0].split()[0])  # current wavelength
                        max_wave = float(data[-1].split()[0])  # maximum wavelength
                        c_nr = 0
                        for d, _ in enumerate(data):
                            cur_wave = float(data[d].split()[0])  # current wavelength
                            if cur_wave < 0.7 * max_wave:
                                break
                            c_nr += 1
                        nl = float(data[c_nr].split()[1])  # lower n value
                        nh = float(data[-1].split()[1])  # higher n value
                        kl = float(data[c_nr].split()[2])  # lower k value
                        kh = float(data[-1].split()[2])  # higher k value

                        # calculate interpolation
                        fac = np.log(wave/max_wave)/np.log(cur_wave/max_wave)
                        ref_index[w, 0] = np.exp(np.log(nh) + fac * np.log(nl/nh))
                        if kl <= 0 or kh < 0:
                            ref_index[w, 1] = 0
                        else:
                            ref_index[w, 1] = np.exp(np.log(kh) + fac * np.log(kl/kh))
                    else:
                        # non-conducting interpolation, linear decreasing k, constant n
                        ref_index[w, 0] = float(data[-1].split()[1])
                        ref_index[w, 1] = float(data[-1].split()[2])*float(data[-1].split()[0])/wave

            # stop searching the files
            break
    else:
        # if no data was found, throw error
        wrt.write_status(
            "WARN",
            'The following cloud particle bulk species is currently ' +
            'not supported or misspelled: ' + species + '. Vacuum values ' +
            'are used instead'
            )

    # return refractive indices
    return ref_index


def patch_mix_opa_tot(self, abundances, mmw, gravity,
                sigma_lnorm=None, fsed=None, Kzz=None,
                radius=None, gray_opacity=None,
                add_cloud_scat_as_abs=None, clouds=None):
    """
    This function is used to update petitRADTRANS to include clouds.
    """
    # Combine total line opacities,
    # according to mass fractions (abundances),
    # also add continuum opacities, i.e. clouds, CIA...

    self.scat = False
    self.mmw = mmw
    for i_spec in range(len(self.line_species)):
        self.line_abundances[:,i_spec] = abundances[self.line_species[i_spec]]
    self.continuum_opa = np.zeros_like(self.continuum_opa)
    self.continuum_opa_scat = np.zeros_like(self.continuum_opa_scat)
    self.continuum_opa_scat_emis = np.zeros_like(self.continuum_opa_scat_emis)

    # Calc. CIA opacity
    if self.H2H2CIA:
        self.continuum_opa = self.continuum_opa + \
                             self.interpolate_cia(self.cia_h2h2_lambda, \
                                                  self.cia_h2h2_temp, self.cia_h2h2_alpha_grid, \
                                                  abundances['H2'],2.)

    if self.H2HeCIA:
        self.continuum_opa = self.continuum_opa + \
                             self.interpolate_cia(self.cia_h2he_lambda, \
                                                  self.cia_h2he_temp, self.cia_h2he_alpha_grid, \
                                                  np.sqrt(abundances['H2']*abundances['He']),np.sqrt(8.))

    if self.N2N2CIA:
        self.continuum_opa = self.continuum_opa + \
                             self.interpolate_cia(self.cia_n2n2_lambda, \
                                                  self.cia_n2n2_temp, self.cia_n2n2_alpha_grid, \
                                                  abundances['N2'],28.)

    if self.O2O2CIA:
        self.continuum_opa = self.continuum_opa + \
                             self.interpolate_cia(self.cia_o2o2_lambda, \
                                                  self.cia_o2o2_temp, self.cia_o2o2_alpha_grid, \
                                                  abundances['O2'],32.)

    if self.N2O2CIA:
        self.continuum_opa = self.continuum_opa + \
                             self.interpolate_cia(self.cia_n2o2_lambda, \
                                                  self.cia_n2o2_temp, self.cia_n2o2_alpha_grid, \
                                                  np.sqrt(abundances['N2']*abundances['O2']),np.sqrt(896.))

    if self.CO2CO2CIA:
        self.continuum_opa = self.continuum_opa + \
                             self.interpolate_cia(self.cia_co2co2_lambda, \
                                                  self.cia_co2co2_temp, self.cia_co2co2_alpha_grid, \
                                                  abundances['CO2'],44.)


    # Calc. H- opacity
    if self.Hminus:
        self.continuum_opa = \
            self.continuum_opa + pyi.hminus_opacity(self.lambda_angstroem, \
                                                    self.border_lambda_angstroem, \
                                                    self.temp, self.press, mmw, abundances)

    # Add mock gray cloud opacity here
    if self.gray_opacity != None:
        self.continuum_opa = self.continuum_opa + self.gray_opacity

    # Add cloud opacity here, will modify self.continuum_opa
    if int(len(self.cloud_species)) > 0:
        self.scat = True
        self.calc_cloud_opacity(abundances, mmw, gravity, \
                                sigma_lnorm, fsed, Kzz, radius, \
                                add_cloud_scat_as_abs)

    # Calculate rayleigh scattering opacities
    if len(self.rayleigh_species) != 0:
        self.scat = True
        self.add_rayleigh(abundances)
    # Add gray cloud deck
    if (self.Pcloud != None):
        self.continuum_opa[:,self.press>self.Pcloud*1e6] += 1e99
    # Add power law opacity
    if (self.kappa_zero != None):
        self.scat = True
        wlen_micron = nc.c/self.freq/1e-4
        scattering_add = self.kappa_zero * \
                         (wlen_micron/0.35)**self.gamma_scat
        add_term = np.repeat(scattering_add[None], \
                             int(len(self.press)), axis = 0).transpose()
        self.continuum_opa_scat += \
            add_term
        if self.do_scat_emis:
            self.continuum_opa_scat_emis += \
                add_term

    # Interpolate line opacities, combine with continuum oacities
    self.line_struc_kappas = fi.mix_opas_ck(self.line_abundances, \
                                            self.line_struc_kappas,self.continuum_opa)

    # Similar to the line-by-line case below, if test_ck_shuffle_comp is
    # True, we will put the total opacity into the first species slot and
    # then carry the remaining radiative transfer steps only over that 0
    # index.
    if (self.mode == 'c-k') and self.test_ck_shuffle_comp:
        #stamps = []
        #import time
        #stamps.append(time.clock())
        #self.combine_opas_shuffle_ck()
        #stamps.append(time.clock())
        '''
        line_struc_kappas_comb = \
          fs.combine_opas_sample_ck(self.line_struc_kappas, \
                                      self.g_gauss, self.w_gauss, \
                                      160)
        '''
        self.line_struc_kappas[:, :, 0, :] = \
            fs.combine_opas_sample_ck(self.line_struc_kappas, \
                                      self.g_gauss, self.w_gauss, \
                                      1000)
        #stamps.append(time.clock())
        #self.combine_opas_shuffle_ck()
        #stamps.append(time.clock())
        #print("Times", np.diff(stamps), \
        #          np.diff(stamps)/np.sum(np.diff(stamps))*100)
        #sys.exit(1)

    # In the line-by-line case we can simply
    # add the opacities of different species
    # in frequency space. All opacities are
    # stored in the first species index slot
    if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
        self.line_struc_kappas[:,:,0,:] = \
            np.sum(self.line_struc_kappas, axis = 2)

    """
    This pit of code is added so clouds can be included. The rest is equivalent to petitRADTRANS
    """
    lsk = self.line_struc_kappas
    self.line_struc_kappas = np.zeros((lsk.shape[0], lsk.shape[1], lsk.shape[2]+1, lsk.shape[3]))
    self.line_struc_kappas[:, :, :-1, :] = lsk
    for j in range(len(self.line_struc_kappas[:, 0, 0, 0])):
        self.line_struc_kappas[j, :, -1, :] = clouds[0]
    self.continuum_opa_scat_emis += clouds[1]


def patch_calc_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl, \
                sigma_lnorm = None, \
                fsed = None, Kzz = None, radius = None, \
                Pcloud = None, \
                kappa_zero = None, \
                gamma_scat = None, \
                contribution = False, haze_factor = None, \
                gray_opacity = None, variable_gravity=True, clouds=None):
    ''' Method to calculate the atmosphere's transmission radius
    (for the transmission spectrum).

        Args:
            temp:
                the atmospheric temperature in K, at each atmospheric layer
                (1-d numpy array, same length as pressure array).
            abunds:
                dictionary of mass fractions for all atmospheric absorbers.
                Dictionary keys are the species names.
                Every mass fraction array
                has same length as pressure array.
            gravity (float):
                Surface gravity in cgs at reference radius and pressure.
            mmw:
                the atmospheric mean molecular weight in amu,
                at each atmospheric layer
                (1-d numpy array, same length as pressure array).
            P0_bar (float):
                Reference pressure P0 in bar where R(P=P0) = R_pl,
                where R_pl is the reference radius (parameter of this
                method), and g(P=P0) = gravity, where gravity is the
                reference gravity (parameter of this method)
            R_pl (float):
                Reference radius R_pl, in cm.
            sigma_lnorm (Optional[float]):
                width of the log-normal cloud particle size distribution
            fsed (Optional[float]):
                cloud settling parameter
            Kzz (Optional):
                the atmospheric eddy diffusion coeffiecient in cgs untis
                (i.e. :math:`\\rm cm^2/s`),
                at each atmospheric layer
                (1-d numpy array, same length as pressure array).
            radius (Optional):
                dictionary of mean particle radii for all cloud species.
                Dictionary keys are the cloud species names.
                Every radius array has same length as pressure array.
            contribution (Optional[bool]):
                If ``True`` the transmission and emission
                contribution function will be
                calculated. Default is ``False``.
            gray_opacity (Optional[float]):
                Gray opacity value, to be added to the opacity at all
                pressures and wavelengths (units :math:`\\rm cm^2/g`)
            Pcloud (Optional[float]):
                Pressure, in bar, where opaque cloud deck is added to the
                absorption opacity.
            kappa_zero (Optional[float]):
                Scarttering opacity at 0.35 micron, in cgs units (cm^2/g).
            gamma_scat (Optional[float]):
                Has to be given if kappa_zero is definded, this is the
                wavelength powerlaw index of the parametrized scattering
                opacity.
            haze_factor (Optional[float]):
                Scalar factor, increasing the gas Rayleigh scattering
                cross-section.
            variable_gravity (Optional[bool]):
                Standard is ``True``. If ``False`` the gravity will be
                constant as a function of pressure, during the transmission
                radius calculation.
            add_cloud_scat_as_abs (Optional[bool]):
                If ``True``, 20 % of the cloud scattering opacity will be
                added to the absorption opacity, introduced to test for the
                effect of neglecting scattering.
    '''

    self.Pcloud = Pcloud
    self.gray_opacity = gray_opacity
    self.interpolate_species_opa(temp)
    self.haze_factor = haze_factor
    self.kappa_zero = kappa_zero
    self.gamma_scat = gamma_scat
    self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius,clouds=clouds)
    self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution,variable_gravity)