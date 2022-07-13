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


class Interface:
    """
    The GCMtools interfacing class which implements common
    interfacing functionalities to different codes.

    Methods
    -------
    _regrid: regrids data to different resolution

    """
    def _regrid(self, data, dlon, dlat):
        """
        Regrid data to a different resolution

        Parameters
        ----------
        data: xr.DataSet or xr.DataArray
            The data that should be regridded
        dlon: float
            the steps in longitutde (longitude resolution)
        dlat: float
            the steps in latititude (latitude resolution)

        Returns
        -------
        data_interp: xr.DataSet or xr.DataArray
            The data that has been regridded
        """
        lonb = np.arange(-180,180+dlon,dlon)
        latb = np.arange(-90, 90+dlat, dlat)

        lons = 0.5 * (lonb[1:] + lonb[:-1])
        lats = 0.5 * (latb[1:] + latb[:-1])

        return data.interp(lon=lons, lat=lats)


class PACInterface(Interface):
    def to_2D_PAC(self, ds):
        """
        Transform the given GCM dataset to a pseudo-2D PAC
        input file.

        Parameters
        ----------
        ds: DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.

        """
        pass


class pRTInterface(Interface):
    """
    Interface with petitRADTRANS
    """

    def __init__(self, pRT):
        """
        Constructs the Interface and links to pRT.

        Parameters
        ----------
        pRT: petitRADTRANS.Radtrans
            A pRT Radtrans object
        """
        self.chemistry = Chemistry()
        self.pRT = pRT

    def set_data(self, ds, regrid_lowres_phase=False):
        """
        Set the data to be used for the interface

        Parameters
        ----------
        ds: DataSet
            A GCMtools-compatible dataset of a 3D climate simulation.
            Should not include a time coordinate!
        regrid_lowres_phase: bool, optional
            regrid the dataset to a 15deg resolution.
            Needed for Phasecurve
        """
        if regrid_lowres_phase:
            ds = self._regrid(data=ds, dlon=15, dlat=15)

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
        return self.chemistry.chem_from_poorman(temp_key=temp_key, CO=CO, FeH=FeH)

    def calc_phase_spectrum(self, mmw, filename=None, **prt_args):
        """
        Calculate the spectrum for a phasecurve.

        Parameters
        ----------
        mmw: float
            Mean molecular weight (in atomic units).
        filename: str
            path at which the output should be stored.
        prt_args:
            All the args that should be parsed to calc_spectra.
            See the docs of prt_phasecurve for more info on the arguments.
            Needed args:
                Rstar, gravity

        Returns
        -------
        spectra: xr.DataArray
            Dataarray containing the spectrum. The spectrum is normed to the stellar spectrum.

        """
        from prt_phasecurve import calc_spectra
        import petitRADTRANS.nat_cst as nc

        abus = self.chemistry.to_pRT(self.pRT.line_species, self.pRT.press/1e6)
        mu = np.cos(abus[c['lon']] * np.pi / 180.) * np.cos(abus[c['lat']] * np.pi / 180.)
        theta_star = np.arccos(mu) * 180 / np.pi
        pRT_abu_keys = list(set(abus.keys()) - {'T', 'nabla_ad', 'MMW'})

        lon, lat = np.meshgrid(abus[c['lon']], abus[c['lat']])

        theta_list, temp_list, abunds_list = [], [], []
        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            theta_list.append(theta_star.sel(lon=lon_i, lat=lat_i).values)
            temp_list.append(abus[c['T']].sel(lon=lon_i, lat=lat_i).values)
            abunds_list.append({sp: abus[sp].sel(lon=lon_i, lat=lat_i).values for sp in pRT_abu_keys})

        gravity_cgs = prt_args.pop('gravity', self.ds.attrs.get(c['g'])*100)
        wlen = nc.c / self.pRT.freq / 1e-4
        stellar_spectrum = self._get_stellar_spec(wlen=wlen, **prt_args)

        spectra_raw = calc_spectra(self.pRT, temp=temp_list, gravity=gravity_cgs, mmw=mmw, abunds=abunds_list,
                                    theta_star=theta_list, **prt_args)

        R_star = prt_args.get('Rstar')
        if R_star is None:
            raise ValueError('pRT needs Rstar, please set Rstar to the stellar radius in cm.')

        R_p = self.ds.attrs.get(c['R_p'])
        if R_p is None:
            raise ValueError(f'pRT needs the planetary radius [in m]. Please use this function only with a GCMT processed dataset.')

        spectra_raw = np.array(spectra_raw) * ((R_p*100)/(R_star))**2
        spectra_raw = spectra_raw/stellar_spectrum[np.newaxis, np.newaxis, :]

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

        ph_c = phase_curve(phases=phases, intensity=spectra, lon=lon, lat=lat)

        return xr.DataArray(
            data=ph_c,
            dims=['phase', 'wlen'],
            coords={'phase': phases,
                    'wlen': spectra.wlen})

    def _get_stellar_spec(self, wlen, Tstar, **kwargs):
        """
        Helperfunction from petitRADTRANS that calculates the stellar spectrum from the phoenix spectrum
        """
        from petitRADTRANS.nat_cst import get_PHOENIX_spec_rad
        if Tstar != None:
            spec, _ = get_PHOENIX_spec_rad(Tstar)
            stellar_intensity = np.interp(wlen*1e-4, spec[:,0], spec[:,1])
            # stellar_intensity = stellar_intensity * (Rstar / semimajoraxis) ** 2

            return stellar_intensity



class Chemistry:
    """
    Chemistry class used to deal with different kinds of chemical models.
    """
    abunds = xr.Dataset()

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
                dims=[c['lon'],c['lat'],c['Z']],
                coords={c['lon']: self.ds[c['lon']],
                        c['lat']: self.ds[c['lat']],
                        c['Z']: self.ds[c['Z']]})})

        for lon in self.ds.lon:
            for lat in self.ds.lat:
                temp_i = self.ds[temp_key].sel(lon=lon, lat=lat)
                abus = interpol_abundances(COs, FeHs, temp_i, p)
                for sp in var_names:
                    self.abunds[sp].loc[{c['lon']:lon, c['lat']:lat}] = abus[sp]
                self.abunds[c['T']].loc[{c['lon']:lon, c['lat']:lat}] = temp_i

        self.abunds.attrs.update({'p_unit':p_unit})

    def to_pRT(self, prt_species, p_prt):
        """
        Format chemistry to match petitRADTRANS requirements.

        Parameters
        ----------
        prt_species: list
            List of species names used in petitRADTRANS
        p_prt: list
            List of species names used in petitRADTRANS

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
            p_prt = p_prt*1e5

        prt_abu = prt_abu.interp(Z=p_prt)

        return prt_abu
