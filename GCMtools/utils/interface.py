# ==============================================================
#                       GCMtools Interface Class
# ==============================================================
#  This class incorporates interfacing methods for GCMtools
#  datasets. The interfacing options are driven by the frequent
#  needs of users to transform GCM output to formats readable by
#  e.g.: - petitRADTRANS            (Molliere+2019)
#        - gCMCRT                   (Lee+2022)
#        - 1D and pseudo-2D PAC     (Agundez+2014)
#        - ...
# ==============================================================

import numpy as np
import xarray as xr

from ..core.const import VARNAMES as c


class _Chemistry:
    """
    Chemistry class used to deal with different kinds of chemical models.
    """

    def __init__(self):
        """
        Constructor for the Chemistry class
        """
        self.abunds = xr.Dataset()
        self.ds = None

    def set_data(self, ds):
        """
        Construct the chemistry class

        Parameters
        ----------
        ds: DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
            Should not have a time coordinate anymore.

        """
        self.ds = ds

    def chem_from_poorman(self, temp_key='T', CO=0.55, FeH=0.0):
        """
        Calculate equilibrium abundancies with poorman from pRT

        Parameters
        ----------
        temp_key: str, optional
            The key to the temperature field used for the abundancies.
            Defaults to T.
        CO: float, optional
            The C/O ratio. Currently only one global value allowed. Defaults to 0.55.
        FeH: float, optional
            The metalicity ratio. Currently only one global value allowed. Defaults to 0.0.
        """
        from petitRADTRANS.poor_mans_nonequ_chem import interpol_abundances

        if self.ds is None:
            raise ValueError("Data is missing. Use interface.set_data() first.")

        if c['time'] in self.ds.dims:
            raise ValueError('Dataset should not have a timedimension. Select the timestamp beforehand.')

        p = self.ds[c['Z']].values
        p_unit = self.ds.attrs.get('p_unit')
        if p_unit == 'Pa':
            p = p / 1e5  # convert to bar
        if p_unit not in ['Pa', 'bar']:
            raise NotImplementedError('can currently only deal with Pa or bar')

        COs = np.ones_like(p) * CO
        FeHs = np.ones_like(p) * FeH

        var_names = interpol_abundances(np.array([0.55]), np.array([0.0]), np.array([100]), np.array([0.1])).keys()
        for key in [c['T'], *var_names]:
            self.abunds.update({key: xr.DataArray(
                data=np.empty((len(self.ds[c['lon']]), len(self.ds[c['lat']]), len(self.ds[c['Z']]))),
                dims=[c['lon'], c['lat'], c['Z']],
                coords={c['lon']: self.ds[c['lon']],
                        c['lat']: self.ds[c['lat']],
                        c['Z']: self.ds[c['Z']]})})

        for lon in self.ds.lon:
            for lat in self.ds.lat:
                temp_i = self.ds[temp_key].sel(lon=lon, lat=lat)
                abus = interpol_abundances(COs, FeHs, temp_i, p)
                for sp in var_names:
                    self.abunds[sp].loc[{c['lon']: lon, c['lat']: lat}] = abus[sp]
                self.abunds[c['T']].loc[{c['lon']: lon, c['lat']: lat}] = temp_i

        self.abunds.attrs.update({'p_unit': p_unit})

    def to_pRT(self, prt_species, p_prt):
        """
        Format chemistry to match petitRADTRANS requirements.

        Parameters
        ----------
        prt_species: list
            List of species names used in petitRADTRANS
        p_prt: list
            Pressure in bar, as used in petitRADTRANS

        Returns
        -------
        prt_abu: Dataset
            abundances ready for prt
        """
        prt_abu = self.abunds.copy()
        for sp_raw in prt_species:
            sp = sp_raw.split('_')[0]
            if sp not in self.abunds:
                raise ValueError(f'We miss chemistry data for {sp}')

            prt_abu.update({sp_raw: self.abunds[sp]})

        p_unit = self.abunds.attrs.get('p_unit')
        if p_unit not in ['Pa', 'bar']:
            raise NotImplementedError('can currently only deal with Pa or bar')

        if p_unit == 'Pa':
            p_prt = p_prt * 1e5

        prt_abu = prt_abu.interp(Z=p_prt)

        return prt_abu


class Interface:
    """
    The GCMtools interfacing class which implements common
    interfacing functionalities to different codes.

    Methods
    -------
    set_data: Function that handles the input from GCMT
    chem_from_poorman: Function that calculates chemistry based on the poorman code from pRT
    """

    def __init__(self, gcmt):
        """
        Constructor for the Interface class

        Parameters
        ----------
        gcmt: GCMTools Object
        """
        self.gcmt = gcmt
        self.chemistry = _Chemistry()
        self.ds = None

    def set_data(self, time, tag=None, regrid_lowres=False):
        """
        Set the data to be used for the interface

        Parameters
        ----------
        time: int
            timestep to be used
        tag: str
            tag of the model to be used
        regrid_lowres: bool, optional
            Can be useful, if your GCMT uses a very detailed grid
        """
        self._set_data_common(time, tag=tag, regrid_lowres=False)

    def _set_data_common(self, time, tag=None, regrid_lowres=False):
        """
        Set the data to be used for the interface

        Parameters
        ----------
        time: int
            timestep to be used
        tag: str
            tag of the model to be used
        regrid_lowres: bool, optional
            Can be useful, if your GCMT uses a very detailed grid
        """
        ds = self.gcmt.get_one_model(tag).sel(time=time)

        if regrid_lowres:
            dlon = 15
            dlat = 15

            lonb = np.arange(-180, 180 + dlon, dlon)
            latb = np.arange(-90, 90 + dlat, dlat)

            lons = 0.5 * (lonb[1:] + lonb[:-1])
            lats = 0.5 * (latb[1:] + latb[:-1])

            ds = ds.interp(lon=lons, lat=lats)

        self.ds = ds
        self.chemistry.set_data(ds)

    def chem_from_poorman(self, temp_key='T', CO=0.55, FeH=0.0):
        """
        Calculate equilibrium abundancies with poorman from pRT

        Parameters
        ----------
        temp_key: str, optional
            The key to the temperature field used for the abundancies.
            Defaults to T.
        CO: float, optional
            The C/O ratio. Currently only one global value allowed. Defaults to 0.55.
        FeH: float, optional
            The metalicity ratio. Currently only one global value allowed. Defaults to 0.0.
        """
        if self.ds is None:
            raise ValueError("Data is missing. Use interface.set_data() first.")

        return self.chemistry.chem_from_poorman(temp_key=temp_key, CO=CO, FeH=FeH)


class pRTInterface(Interface):
    """
    Interface with petitRADTRANS
    """

    def __init__(self, gcmt, pRT):
        """
        Constructs the Interface and links to pRT.

        Parameters
        ----------
        gcmt: GCMTools Object
            A GCMTools object that is linked to the interface
        pRT: petitRADTRANS.Radtrans
            A pRT Radtrans object
        """
        super().__init__(gcmt)
        self.pRT = pRT

    def set_data(self, time, tag=None, regrid_lowres=False):
        """
        Set the data to be used for the interface

        Parameters
        ----------
        time: int
            timestep to be used
        tag: str
            tag of the model to be used
        regrid_lowres: bool, optional
            Can be useful, if your GCMT uses a very detailed grid
        """
        self._set_data_common(time, tag=tag, regrid_lowres=False)

        if self.ds.p_unit == 'bar':
            press = self.ds.Z.values
        elif self.ds.p_unit == 'Pa':
            press = self.ds.Z.values/1e5
        else:
            raise NotImplementedError('only pressure units in Pa and bar are implemented at the moment')

        # be careful here: petitRADTRANS operates top->bot in bar
        self.pRT.setup_opa_structure(np.sort(press))

    def calc_phase_spectrum(self, mmw, Rstar, Tstar, semimajoraxis, gravity=None, filename=None, normalize=True,
                            **prt_args):
        """
        Calculate the spectrum for a phasecurve.

        Parameters
        ----------
        mmw: float or 1D-array
            Mean molecular weight (in atomic units). Will be globally uniform if float or horizonatally uniform if 1D.
        Rstar: float
            stellar radius in !cm!.
        Tstar: float
            Temperature of the hoststar.
        semimajoraxis: float
            The distance between hoststar and planet in !cm!
        gravity: float
            surface gravity in !cgs!. Will default to the value provided by GCMT.
        filename: str
            path at which the output should be stored.
        normalize: bool
            Decide if you want the spectrum to be normalized to the star. The code will then:
            1. correct intensity for the ratio between solar radius and planetary radius
            2. devide by stellar spectrum
        prt_args:
            All the args that should be parsed to calc_spectra.
            See the docs of prt_phasecurve for more info on the arguments.

        Returns
        -------
        spectra: xr.DataArray
            Dataarray containing the spectrum. The spectrum is normed to the stellar spectrum.

        """
        from prt_phasecurve import calc_spectra
        import petitRADTRANS.nat_cst as nc

        abus = self.chemistry.to_pRT(self.pRT.line_species, self.pRT.press / 1e6)
        mu = np.cos(abus[c['lon']] * np.pi / 180.) * np.cos(abus[c['lat']] * np.pi / 180.)
        theta_star = np.arccos(mu) * 180 / np.pi
        pRT_abu_keys = list(set(abus.keys()) - {c['T'], 'nabla_ad', 'MMW'})

        lon, lat = np.meshgrid(abus[c['lon']], abus[c['lat']])

        if (lon.shape[0] * lon.shape[1]) > 300:
            print(
                'WARNING: Calculating a phasecurve on a fine grid takes very long. A resolution of 15 degrees is usually sufficient.')

        theta_list, temp_list, abunds_list = [], [], []
        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            theta_list.append(theta_star.sel(lon=lon_i, lat=lat_i).values)
            temp_list.append(abus[c['T']].sel(lon=lon_i, lat=lat_i).values)
            abunds_list.append({sp: abus[sp].sel(lon=lon_i, lat=lat_i).values for sp in pRT_abu_keys})

        if gravity is None:
            gravity = self.ds.attrs.get(c['g']) * 100

        wlen = nc.c / self.pRT.freq / 1e-4
        stellar_spectrum = self._get_stellar_spec(wlen=wlen, Tstar=Tstar)
        mmw = np.ones_like(self.pRT.press) * mmw  # broadcast if needed

        spectra_raw = calc_spectra(self.pRT, temp=temp_list, gravity=gravity, mmw=mmw, abunds=abunds_list,
                                   theta_star=theta_list, Tstar=Tstar, Rstar=Rstar, semimajoraxis=semimajoraxis,
                                   **prt_args)
        spectra_raw = np.array(spectra_raw)

        R_p = self.ds.attrs.get(c['R_p'])
        if R_p is None:
            raise ValueError(
                f'pRT needs the planetary radius [in m]. Please use this function only with a GCMT processed dataset.')

        if normalize:
            # correct by (R_p/R_s)**2 and norm to stellar spectrum:
            spectra_raw = spectra_raw * ((R_p * 100) / (Rstar)) ** 2
            spectra_raw = spectra_raw / stellar_spectrum[np.newaxis, np.newaxis, :]

        Nmus = spectra_raw.shape[1]
        spectra = xr.DataArray(
            data=np.empty((len(abus[c['lon']]), len(abus[c['lat']]), Nmus, len(wlen))),
            dims=[c['lon'], c['lat'], 'imu', 'wlen'],
            coords={c['lon']: abus[c['lon']],
                    c['lat']: abus[c['lat']],
                    'imu': range(Nmus),
                    'wlen': wlen})

        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            spectra.loc[{c['lon']: lon_i, c['lat']: lat_i}] = spectra_raw[i]

        if filename is not None:
            spectra.to_netcdf(filename)

        return spectra

    def phase_curve(self, phases, spectra=None, filename=None):
        """
        Do the diskintegration of the spectrum to yield the phasecurve

        Parameters
        ----------
        phases (array(P)):
            List of phases at which the phasecurve should be evaluated
        spectra: str
            The spectrum generated by calc_phase_spectrum
        filename: str
            Alternatively give a filename where the spectrum is saved

        Returns
        -------
        phasecurve: xr.DataArray
            DataArray containing the phasecurveinformation

        """
        from prt_phasecurve import phase_curve

        if spectra is None and filename is None:
            raise ValueError('please provide a file with the spectrum or a spectrum produced by calc_phase_spectrum')

        if filename is not None and spectra is None:
            spectra = xr.open_dataarray(filename)

        lon, lat = np.meshgrid(spectra[c['lon']], spectra[c['lat']])

        intensity, lat_list, lon_list = [], [], []
        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            lat_list.append(lat_i)
            lon_list.append(lon_i)
            intensity.append(spectra.sel(lon=lon_i, lat=lat_i).values)

        ph_c = phase_curve(phases=phases, intensity=intensity, lon=lon_list, lat=lat_list)

        return xr.DataArray(
            data=ph_c,
            dims=['phase', 'wlen'],
            coords={'phase': phases,
                    'wlen': spectra.wlen})

    def _get_stellar_spec(self, wlen, Tstar):
        """
        Helperfunction from petitRADTRANS that calculates the stellar spectrum from the phoenix spectrum
        """
        from petitRADTRANS.nat_cst import get_PHOENIX_spec_rad
        if Tstar != None:
            spec, _ = get_PHOENIX_spec_rad(Tstar)
            stellar_intensity = np.interp(wlen * 1e-4, spec[:, 0], spec[:, 1])
            return stellar_intensity

# class PACInterface(Interface):
#     def to_2D_PAC(self, ds):
#         """
#         Transform the given GCM dataset to a pseudo-2D PAC
#         input file.
#
#         Parameters
#         ----------
#         ds: DataSet
#             A GCMtools-compatible dataset of a 3D climate simulation.
#
#         """
#         pass
