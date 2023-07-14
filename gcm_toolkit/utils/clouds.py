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
from glob import glob
import numpy as np
from scipy.optimize import minimize
from ..core import writer as wrt


def cloud_opacities(wavelengths, cloud_radius, cloud_abundances, cloud_particle_density,
                    volume_fraction, use_bruggemann=False):
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
    use_bruggemann : bool
        If this flag is set to true, bruggemann mixing is used. This is
        more accurate but takes much longer to calculate.


    Returns
    -------
    Qabs : np.ndarray[h, w]
        Absorption coefficients
    Qsca : np.ndarray[h, w]
        Scattering coefficients
    c_sec : np.ndarray[h]
        cross sections
    """

    # load mie calculation library
    import miepython

    # check input
    if not len(cloud_radius) == len(cloud_abundances) == len(cloud_particle_density):
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
    eff = eff_func(volume_fraction, wavelengths, use_bruggemann=use_bruggemann)

    # loop over all inputs:
    for wav in range(len_w):
        for hig in range(len_h):
            # prepare input parameters for mie calculation
            m_eff = complex(eff[wav, hig, 0], -eff[wav, hig, 1])
            x_fac = 2 * np.pi * cloud_radius[hig] / wavelengths[wav]

            # calculate mie theory
            qe_temp, qs_temp, _, _ = miepython.mie(m_eff, x_fac)

            # save data
            qabs[hig, wav] = qe_temp - qs_temp
            qsca[hig, wav] = qs_temp

    # calculate cloud particle cross section
    c_sec = cloud_abundances * 3 / 4 / (cloud_radius * cloud_particle_density)

    # all done
    return qabs, qsca, c_sec


def eff_func(volume_fraction, wavelength, use_bruggemann=False):
    """
    Get the effective real and imaginary refractory index of a mixed
    cloud particle grain. As default, the LLL method is used.

    Parameters
    ----------
    volume_fraction: dict[Species: np.ndarray[h]]
        Volume fractions of cloud species. All species included must be
        supported by the cloud_nk_data function.
    wavelength: np.ndarray[w]
        The wavelength at which the refractory index should be given [m].
    use_bruggemann : bool
        If this flag is set to true, bruggemann mixing is used. This is
        more accurate but takes much longer to calculate.

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
        for hig in range(len_h):
            work[k, hig, :, 1:] = temp

    # get all volume fractions
    for wav in range(len_w):
        for hig in range(len_h):
            for k, key in enumerate(volume_fraction):
                work[k, hig, wav, 0] = volume_fraction[key][hig]

    # loop over all height points
    for wav in range(len_w):
        for hig in range(len_h):
            if use_bruggemann:
                wrt.write_status('INFO', 'Bruggemann progress: ' + str(round((hig + wav*len_h)/len_h/len_w*100, 1)) + '  %', end='\r')
                # initial guess using linear approximation
                eff_0 = np.zeros((2,))
                eff_0[0] = sum(work[:, hig, wav, 0] * work[:, hig, wav, 1])
                eff_0[1] = sum(work[:, hig, wav, 0] * work[:, hig, wav, 2])

                # calculate effective medium theory with Bruggemann minimization
                eff[wav, hig] = minimize(_func, eff_0, args=work[:, hig, wav, :]).x
            else:
                # LLL method
                val = complex(0, 0)
                for k, key in enumerate(volume_fraction):
                    val += work[k, hig, wav, 0] * complex(work[k, hig, wav, 1],
                                                          work[k, hig, wav, 2]) ** (1 / 3)
                val = val ** 3
                tmp = np.sqrt(val.real ** 2 + val.imag ** 2)
                eff[wav, hig, 0] = np.sqrt(0.5 * (val.real + tmp))
                eff[wav, hig, 1] = np.sqrt(0.5 * (-val.real + tmp)) * val.imag / np.abs(val.imag)

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
        ava = work[i, 1] ** 2 - work[i, 2] ** 2
        bva = 2 * work[i, 1] * work[i, 2]
        cva = eff[0]
        dva = eff[1]

        # calcualte vol_frac*(e - eff)/(e + 2*eff) for complex numbers
        temp = ava ** 2 + 4 * ava * cva + bva ** 2 + 4 * bva * dva + 4 * (cva ** 2 + dva ** 2)
        val[0] += work[i, 0] * (ava ** 2 + ava * cva + bva ** 2 + bva * dva - 2 * (cva ** 2 + dva ** 2)) / temp
        val[1] += work[i, 0] * (3 * bva * cva - 3 * ava * dva) / temp

    return val[0] ** 2 + val[1] ** 2


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
    wave_temp = wavelength * 1e6

    # prepare output
    if not isinstance(wave_temp, np.ndarray):
        wave_temp = np.asarray([wave_temp])
    ref_index = np.zeros((len(wave_temp), 2))

    # find species in files
    for fil in files:
        if species in fil:

            # get data
            content = open(fil, 'r').readlines()
            header = content[0].split()[1]
            data = content[5:]

            # loop over all wavelengths:
            for wav, wave in enumerate(wave_temp):
                # if desired wavelength is smaller than data, use the smallest
                # wavelength data available.
                if wave < float(data[0].split()[0]):
                    ref_index[wav, 0] = float(data[0].split()[1])
                    ref_index[wav, 1] = float(data[0].split()[2])
                    continue

                # if wavelength is within range, log-log interpolation
                for dnr, _ in enumerate(data):
                    cur_wave = float(data[dnr].split()[0])  # current wavelength
                    if wave < cur_wave:
                        nlo = float(data[dnr - 1].split()[1])  # lower n value
                        nhi = float(data[dnr].split()[1])  # higher n value
                        klo = float(data[dnr - 1].split()[2])  # lower k value
                        khi = float(data[dnr].split()[2])  # higher k value
                        prev_wave = float(data[dnr - 1].split()[0])  # previous wavelength

                        # calculate interpolation
                        fac = np.log(wave / prev_wave) / np.log(cur_wave / prev_wave)
                        ref_index[wav, 0] = np.exp(np.log(nlo) + fac * np.log(nhi / nlo))
                        if klo <= 0 or khi < 0:
                            ref_index[wav, 1] = 0
                        else:
                            ref_index[wav, 1] = np.exp(np.log(klo) + fac * np.log(khi / klo))
                        break
                else:
                    # if wavelength is outside of range, extrapolate
                    if header == '.True.':
                        # conducting interpolation, log-log interpolation
                        # find 70% wavelength
                        cur_wave = float(data[0].split()[0])  # current wavelength
                        max_wave = float(data[-1].split()[0])  # maximum wavelength
                        c_nr = 0
                        for dnr, _ in enumerate(data):
                            cur_wave = float(data[dnr].split()[0])  # current wavelength
                            if cur_wave < 0.7 * max_wave:
                                break
                            c_nr += 1
                        nlo = float(data[c_nr].split()[1])  # lower n value
                        nhi = float(data[-1].split()[1])  # higher n value
                        klo = float(data[c_nr].split()[2])  # lower k value
                        khi = float(data[-1].split()[2])  # higher k value

                        # calculate interpolation
                        fac = np.log(wave / max_wave) / np.log(cur_wave / max_wave)
                        ref_index[wav, 0] = np.exp(np.log(nhi) + fac * np.log(nlo / nhi))
                        if klo <= 0 or khi < 0:
                            ref_index[wav, 1] = 0
                        else:
                            ref_index[wav, 1] = np.exp(np.log(khi) + fac * np.log(klo / khi))
                    else:
                        # non-conducting interpolation, linear decreasing k, constant n
                        ref_index[wav, 0] = float(data[-1].split()[1])
                        ref_index[wav, 1] = float(data[-1].split()[2]) * \
                                            float(data[-1].split()[0]) / wave

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


def patch_cloud_mix_opa(self, clouds=None):
    """
    This function is used to update petitRADTRANS to include clouds.
    """
    if clouds is not None:
        lsk = self.line_struc_kappas
        self.line_struc_kappas = np.zeros((lsk.shape[0], lsk.shape[1],
                                           lsk.shape[2] + 1, lsk.shape[3]))
        self.line_struc_kappas[:, :, :-1, :] = lsk
        self.line_struc_kappas[:, :, -1, :] = np.ones_like(self.line_struc_kappas[:, :, -1, :]) \
                                              * clouds[0][np.newaxis, :, :]
        self.continuum_opa_scat_emis += clouds[1]


def patch_delete_clouds(self, clouds=None):
    """
    This function is used to update petitRADTRANS to include clouds.
    """
    if clouds is not None:
        self.line_struc_kappas = self.line_struc_kappas[:, :, :-1, :]


def patch_calc_transm(self, temp, abunds, gravity, mmw, P0_bar, R_pl,
                      sigma_lnorm=None,
                      fsed=None, Kzz=None, radius=None,
                      Pcloud=None,
                      kappa_zero=None,
                      gamma_scat=None,
                      contribution=False, haze_factor=None,
                      gray_opacity=None, variable_gravity=True, clouds=None):
    """
    This function is used to update petitRADTRANS to include clouds.
    """

    self.Pcloud = Pcloud
    self.gray_opacity = gray_opacity
    self.interpolate_species_opa(temp)
    self.haze_factor = haze_factor
    self.kappa_zero = kappa_zero
    self.gamma_scat = gamma_scat
    self.mix_opa_tot(abunds, mmw, gravity, sigma_lnorm, fsed, Kzz, radius)
    self.clouds_mix_opa(clouds)
    self.calc_tr_rad(P0_bar, R_pl, gravity, mmw, contribution, variable_gravity)
    self.delete_clouds(clouds)
