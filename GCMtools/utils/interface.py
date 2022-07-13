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
            See the docs of prt_phasecurve for more info on the arguments
        """
        from prt_phasecurve import calc_spectra

        mu = np.cos(self.ds[c['lon']] * np.pi / 180.) * np.cos(self.ds[c['lat']] * np.pi / 180.)
        theta_star = np.arccos(mu) * 180 / np.pi
        abus = self.chemistry.abunds.to_pRT(self.pRT.line_species, self.pRT.press/1e6)

        lon, lat = np.meshgrid(self.ds[c['lon']], self.ds[c['lat']])

        theta_list, temp_list, abunds_list = [], [], []
        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            theta_list.append(np.array(theta_star.sel(lon=lon_i, lat=lat_i)))
            temp_list.append(np.array(abus[c['T']].sel(lon=lon_i, lat=lat_i)))
            abunds_list.append(abus.sel(lon=lon_i, lat=lat_i).to_dict())

        spectra_list = calc_spectra(self.pRT, temp=temp_list, mmw=mmw, abunds=abunds_list, theta_star=theta_list, **prt_args)

        spectra = xr.DataArray(
            data=np.empty((len(abus[c['lon']]), len(abus[c['lat']]), len(self.pRT.wlen))),
            dims=[c['lon'], c['lat'], 'wlen'],
            coords={c['lon']: abus[c['lon']],
                    c['lat']: abus[c['lat']],
                    'wlen': self.pRT.wlen})

        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            spectra.loc[{c['lon']: lon_i, c['lat']: lat_i}] = spectra_list[i]

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
        """
        from prt_phasecurve import phase_curve

        if spectra is None and filename is None:
            raise ValueError('please provide a file with the spectrum or a spectrum produced by calc_phase_spectrum')

        if filename is not None and spectra is None:
            spectra = xr.open_dataarray(filename)

        return phase_curve(phases=phases, intensity=spectra, lon=self.ds[c['lat']], lat=self.ds[c['lat']])


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

        var_names = interpol_abundances([0.55], [0.0], [100], [1]).keys()
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

            prt_abu.update({sp_raw: sp})

        p_unit = self.abunds.attrs.get('p_unit')
        if p_unit not in ['Pa', 'bar']:
            raise NotImplementedError('can currently only deal with Pa or bar')

        if p_unit == 'Pa':
            p_prt = p_prt*1e5

        prt_abu = prt_abu.interp(Z=p_prt)

        return prt_abu
