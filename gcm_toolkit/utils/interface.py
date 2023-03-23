"""
==============================================================
                      gcm_toolkit Interface Class
==============================================================
 This class incorporates interfacing methods for gcm_toolkit
 datasets. The interfacing options are driven by the frequent
 needs of users to transform GCM output to formats readable by
 e.g.: - petitRADTRANS            (Molliere+2019)
       - gCMCRT                   (Lee+2022)
       - 1D and pseudo-2D PAC     (Agundez+2014)
       - ...
==============================================================
"""
import math
import types
import numpy as np
import xarray as xr

from ..core import writer as wrt
from ..core.const import VARNAMES as c
from ..utils.passport import is_the_data_cloudy
from ..utils.clouds import cloud_opacities, patch_cloud_mix_opa, patch_calc_transm, \
                           patch_delete_clouds


class _Chemistry:
    """
    Chemistry class used to deal with different kinds of chemical models.
    """

    def __init__(self):
        """
        Constructor for the Chemistry class
        """
        self.abunds = xr.Dataset()
        self.dsi = None

    def set_data(self, dsi):
        """
        Construct the chemistry class

        Parameters
        ----------
        dsi: DataSet
            A gcm_toolkit-compatible dataset of a 3D climate simulation.
            Should not have a time coordinate anymore.

        """
        self.dsi = dsi

    def chem_from_poorman(self, temp_key="T", co_ratio=0.55, feh_ratio=0.0):
        """
        Calculate equilibrium abundancies with poorman from pRT

        Parameters
        ----------
        temp_key: str, optional
            The key to the temperature field used for the abundancies.
            Defaults to T.
        co_ratio: float, optional
            The C/O ratio. Currently only one global value allowed. Defaults to 0.55.
        feh_ratio: float, optional
            The metalicity ratio. Currently only one global value allowed.
            Defaults to 0.0.
        """
        from petitRADTRANS.poor_mans_nonequ_chem import interpol_abundances

        if self.dsi is None:
            raise ValueError(
                "Data is missing. Use interface.set_data() first."
            )

        if c["time"] in self.dsi.dims:
            raise ValueError(
                "Dataset should not have a timedimension. "
                + "Select the timestamp beforehand."
            )

        pres = self.dsi[c["Z"]].values
        p_unit = self.dsi.attrs.get("p_unit")
        if p_unit == "Pa":
            pres = pres / 1e5  # convert to bar
        if p_unit not in ["Pa", "bar"]:
            raise NotImplementedError("can currently only deal with Pa or bar")

        co_ratios = np.ones_like(pres) * co_ratio
        feh_ratios = np.ones_like(pres) * feh_ratio

        var_names = interpol_abundances(
            np.array([0.55]), np.array([0.0]), np.array([100]), np.array([0.1])
        ).keys()
        for key in [c["T"], *var_names]:
            self.abunds.update(
                {
                    key: xr.DataArray(
                        data=np.empty(
                            (
                                len(self.dsi[c["lon"]]),
                                len(self.dsi[c["lat"]]),
                                len(self.dsi[c["Z"]]),
                            )
                        ),
                        dims=[c["lon"], c["lat"], c["Z"]],
                        coords={
                            c["lon"]: self.dsi[c["lon"]],
                            c["lat"]: self.dsi[c["lat"]],
                            c["Z"]: self.dsi[c["Z"]],
                        },
                    )
                }
            )

        for lon in self.dsi.lon:
            for lat in self.dsi.lat:
                temp_i = self.dsi[temp_key].sel(lon=lon, lat=lat)
                abus = interpol_abundances(co_ratios, feh_ratios, temp_i, pres)
                for spi in var_names:
                    self.abunds[spi].loc[
                        {c["lon"]: lon, c["lat"]: lat}
                    ] = abus[spi]
                self.abunds[c["T"]].loc[
                    {c["lon"]: lon, c["lat"]: lat}
                ] = temp_i

        self.abunds.attrs.update({"p_unit": p_unit})

    def to_prt(self, prt_species, p_prt):
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
            spi = sp_raw.split("_")[0]
            if spi not in self.abunds:
                raise ValueError(f"We miss chemistry data for {spi}")

            prt_abu.update({sp_raw: self.abunds[spi]})

        p_unit = self.abunds.attrs.get("p_unit")
        if p_unit not in ["Pa", "bar"]:
            raise NotImplementedError("can currently only deal with Pa or bar")

        if p_unit == "Pa":
            p_prt = p_prt * 1e5

        prt_abu = prt_abu.interp(Z=p_prt)

        return prt_abu


class Interface:
    """
    The gcm_toolkit interfacing class which implements common
    interfacing functionalities to different codes.

    Methods
    -------
    set_data: Function that handles the input from GCMT
    chem_from_poorman:
        Function that calculates chemistry based on the
        poorman code from pRT
    """

    def __init__(self, tools):
        """
        Constructor for the Interface class

        Parameters
        ----------
        tools: GCMTools Object
        """
        self.tools = tools
        self.chemistry = _Chemistry()
        self.dsi = None

    def set_data(self, time, tag=None, regrid_lowres=False,
                 terminator_avg=False, lat_points=1, lon_resolution=10):
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
        terminator_avg: bool, optional
            terminator avaraging. This requires lat_points and lon_resolution.
        lat_points: int, optional
            number of equally spaced latitude grid points for each terminator (morning and evening)
        lon_resolution: float, optional
            longitudinal opening angle for terminator avareging
        """

        self._set_data_common(time, tag=tag, regrid_lowres=regrid_lowres,
                              terminator_avg=terminator_avg, lat_points=lat_points,
                              lon_resolution=lon_resolution)

    def _set_data_common(self, time, tag=None, regrid_lowres=False,
                         terminator_avg=False, lat_points=1, lon_resolution=10):
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
        terminator_avg: bool, optional
            terminator avaraging. This requires lat_points and lon_resolution.
        lat_points: int, optional
            number of equally spaced latitude grid points for each terminator (morning and evening)
        lon_resolution: float, optional
            longitudinal opening angle for terminator avareging
        """
        # check input legality
        if sum([regrid_lowres, terminator_avg]) > 1:
            raise ValueError('In the function Interface.set_data, multiple regriding options '
                             'were selected, but only one can be given.')

        dsi = self.tools.get_one_model(tag).sel(time=time)

        if terminator_avg:
            # avarage over longitudinal opening angle. Note, this step also corrects
            # for smaller areas at polar regions
            ds_m = dsi.where((dsi[c['lon']] > -90 - lon_resolution/2) *
                             (dsi[c['lon']] < -90 + lon_resolution/2)).mean(c['lon'])
            ds_e = dsi.where((dsi[c['lon']] > 90 - lon_resolution/2) *
                             (dsi[c['lon']] < 90 + lon_resolution/2)).mean(c['lon'])

            # set latitude step size
            lat_step = 180 / lat_points
            lat_x = np.linspace(-90+lat_step/2, 90-lat_step/2, lat_points)

            # generate new dataset
            ds_transit = xr.Dataset(
                data_vars={},
                coords={
                    'lat': ([c["lat"]], lat_x),
                    'lon': ([c["lon"]], [-90, 90]),
                    'Z_l': ([c["Z_l"]], dsi[c['Z_l']].values),
                    'Z': ([c["Z"]], dsi[c['Z']].values)
                },
                attrs=dsi.attrs
            )

            # avarage in latitude space
            for key in ds_e.keys():
                if key in ['T', 'ClAb', 'ClDs', 'ClDr'] or 'ClVf' in key:
                    # empty array to fill with data
                    tmp = np.zeros((lat_points, 2, len(dsi[c["Z"]].values)))
                    for lat in range(lat_points):
                        tmp[lat, 0, :] = ds_m.where((ds_m[c['lat']] > lat*lat_step - 90) *
                                                    (ds_m[c['lat']] < (lat+1) * lat_step - 90)
                                                    ).mean(c['lat'])[key].values
                        tmp[lat, 1, :] = ds_e.where((ds_e[c['lat']] > lat*lat_step - 90) *
                                                    (ds_e[c['lat']] < (lat+1)*lat_step - 90)
                                                    ).mean(c['lat'])[key].values

                    # saving results
                    ds_transit[key] = ((c['lat'], c['lon'], c['Z_l']), tmp)

            # add mark that data is ready for tranist callcuations
            ds_transit.attrs['transit'] = True

            # return prepared dataset
            dsi = ds_transit

        if regrid_lowres:
            dlon = 15
            dlat = 15

            lonb = np.arange(-180, 180 + dlon, dlon)
            latb = np.arange(-90, 90 + dlat, dlat)

            lons = 0.5 * (lonb[1:] + lonb[:-1])
            lats = 0.5 * (latb[1:] + latb[:-1])

            dsi = dsi.interp(lon=lons, lat=lats)

        self.dsi = dsi
        self.chemistry.set_data(dsi)
        self.chemistry.abunds = xr.Dataset()

    def chem_from_poorman(self, temp_key="T", co_ratio=0.55, feh_ratio=0.0):
        """
        Calculate equilibrium abundancies with poorman from pRT

        Parameters
        ----------
        temp_key: str, optional
            The key to the temperature field used for the abundancies.
            Defaults to T.
        co_ratio: float, optional
            The C/O ratio. Currently only one global value allowed.
            Defaults to 0.55.
        feh_ratio: float, optional
            The metalicity ratio. Currently only one global value allowed.
             Defaults to 0.0.
        """
        if self.dsi is None:
            raise ValueError(
                "Data is missing. Use interface.set_data() first."
            )

        return self.chemistry.chem_from_poorman(
            temp_key=temp_key, co_ratio=co_ratio, feh_ratio=feh_ratio
        )


class PrtInterface(Interface):
    """
    Interface with petitRADTRANS
    """

    def __init__(self, tools, prt):
        """
        Constructs the Interface and links to pRT.

        Parameters
        ----------
        tools: GCMTools Object
            A GCMTools object that is linked to the interface
        pRT: petitRADTRANS.Radtrans
            A pRT Radtrans object
        """
        super().__init__(tools)
        self.prt = prt

    def set_data(self, time, tag=None, regrid_lowres=False,
                 terminator_avg=False, lat_points=1, lon_resolution=10):
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
        terminator_avg: bool, optional
            terminator avaraging. This requires lat_points and lon_resolution.
        lat_points: int, optional
            number of equally spaced latitude grid points for each terminator (morning and evening)
        lon_resolution: float, optional
            longitudinal opening angle for terminator avareging
        """

        self._set_data_common(time, tag=tag, regrid_lowres=regrid_lowres,
                              terminator_avg=terminator_avg, lat_points=lat_points,
                              lon_resolution=lon_resolution)

        if self.dsi.p_unit == "bar":
            press = self.dsi.Z.values
        elif self.dsi.p_unit == "Pa":
            press = self.dsi.Z.values / 1e5
        else:
            raise NotImplementedError(
                "only pressure units in Pa and bar are "
                + "implemented at the moment"
            )

        # be careful here: petitRADTRANS operates top->bot in bar
        self.prt.setup_opa_structure(np.sort(press))

    # pylint: disable=C0103
    def calc_phase_spectrum(
        self,
        mmw,
        Rstar,
        Tstar,
        semimajoraxis,
        gravity=None,
        filename=None,
        normalize=True,
        **prt_args,
    ):
        """
        Calculate the spectrum for a phasecurve.

        Parameters
        ----------
        mmw: float or 1D-array
            Mean molecular weight (in atomic units). Will be globally uniform if
            float or horizonatally uniform if 1D.
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
            Decide if you want the spectrum to be normalized to the star.
            The code will then:
            1. correct intensity for the ratio between solar radius and
                planetary radius
            2. devide by stellar spectrum
        prt_args:
            All the args that should be parsed to calc_spectra.
            See the docs of prt_phasecurve for more info on the arguments.

        Returns
        -------
        spectra: xr.DataArray
            Dataarray containing the spectrum.
            The spectrum is normed to the stellar spectrum.

        """
        from prt_phasecurve import calc_spectra
        import petitRADTRANS.nat_cst as nc

        abus = self.chemistry.to_prt(
            self.prt.line_species, self.prt.press / 1e6
        )
        mui = np.cos(abus[c["lon"]] * np.pi / 180.0) * np.cos(
            abus[c["lat"]] * np.pi / 180.0
        )
        theta_star = np.arccos(mui) * 180 / np.pi
        prt_abu_keys = list(set(abus.keys()) - {c["T"], "nabla_ad", "MMW"})

        lon, lat = np.meshgrid(abus[c["lon"]], abus[c["lat"]])

        if (lon.shape[0] * lon.shape[1]) > 300:
            raise ValueError(
                "WARNING: Calculating a phasecurve on a fine grid takes very"
                " long. "
                + "A resolution of 15 degrees is usually sufficient."
            )

        theta_list, temp_list, abunds_list = [], [], []
        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            theta_list.append(theta_star.sel(lon=lon_i, lat=lat_i).values)
            temp_list.append(abus[c["T"]].sel(lon=lon_i, lat=lat_i).values)
            abunds_list.append(
                {
                    sp: abus[sp].sel(lon=lon_i, lat=lat_i).values
                    for sp in prt_abu_keys
                }
            )

        if gravity is None:
            gravity = self.dsi.attrs.get(c["g"]) * 100

        wlen = nc.c / self.prt.freq / 1e-4
        stellar_spectrum = self._get_stellar_spec(wlen=wlen, t_star=Tstar)
        mmw = np.ones_like(self.prt.press) * mmw  # broadcast if needed

        spectra_raw = calc_spectra(
            self.prt,
            temp=temp_list,
            gravity=gravity,
            mmw=mmw,
            abunds=abunds_list,
            theta_star=theta_list,
            Tstar=Tstar,
            Rstar=Rstar,
            semimajoraxis=semimajoraxis,
            **prt_args,
        )
        spectra_raw = np.array(spectra_raw)

        r_p = self.dsi.attrs.get(c["R_p"])
        if r_p is None:
            raise ValueError(
                "pRT needs the planetary radius [in m]. Please use this"
                " function "
                + "only with a GCMT processed dataset."
            )

        if normalize:
            # correct by (r_p/R_s)**2 and norm to stellar spectrum:
            spectra_raw = spectra_raw * ((r_p * 100) / (Rstar)) ** 2
            spectra_raw = (
                spectra_raw / stellar_spectrum[np.newaxis, np.newaxis, :]
            )

        nmus = spectra_raw.shape[1]
        spectra = xr.DataArray(
            data=np.empty(
                (len(abus[c["lon"]]), len(abus[c["lat"]]), nmus, len(wlen))
            ),
            dims=[c["lon"], c["lat"], "imu", "wlen"],
            coords={
                c["lon"]: abus[c["lon"]],
                c["lat"]: abus[c["lat"]],
                "imu": range(nmus),
                "wlen": wlen,
            },
        )

        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            spectra.loc[{c["lon"]: lon_i, c["lat"]: lat_i}] = spectra_raw[i]

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
            raise ValueError(
                "please provide a file with the spectrum or a spectrum"
                " produced by calc_phase_spectrum"
            )

        if filename is not None and spectra is None:
            spectra = xr.open_dataarray(filename)

        lon, lat = np.meshgrid(spectra[c["lon"]], spectra[c["lat"]])

        intensity, lat_list, lon_list = [], [], []
        for i, lon_i in enumerate(np.array(lon).flat):
            lat_i = np.array(lat).flat[i]
            lat_list.append(lat_i)
            lon_list.append(lon_i)
            intensity.append(spectra.sel(lon=lon_i, lat=lat_i).values)

        ph_c = phase_curve(
            phases=phases, intensity=intensity, lon=lon_list, lat=lat_list
        )

        return xr.DataArray(
            data=ph_c,
            dims=["phase", "wlen"],
            coords={"phase": phases, "wlen": spectra.wlen},
        )

    def calc_transit_spectrum(
            self,
            mmw,
            gravity=None,
            rplanet=None,
            pressure_0=None,
            mass_frac=None,
            clouds=None,
            use_bruggemann=False
    ):
        """
        Calculate the transit spectrum. This function avarages T-p profiles in
        the terminator region and uses poorman equilbriums chemistry.

        Parameters
        ----------
        mmw: float or 1D-array
            Mean molecular weight (in atomic units). Will be globally uniform if
            float or horizonatally uniform if 1D.
        gravity: float
            surface gravity in !cgs!. Will default to the value provided by GCMT.
        rplanet: float
            planet radius in !cm!. Will default to the value provided by GCMT.
        pressure_0: float
            pressure level of the gravity and radius. Will default to the value provided by GCMT.
        mass_frac: List[List[Dict]], optional
            The mass fractions for each t_p point. Structure is as follows: mass_frac[i][j]['key']
            - i: [2 elements] 0 for morning terminator, 1 for evening terminator
            - j: [lat_points elements] latitude point starting from lowest
                 (closeset to lat = -90) to highest (closest to lat = 90)
            - dict must contain all mass fractions  elements given as opacity species to prt
        clouds: Union[np.ndarray[i, j, h, w, k], bool], optional
            The cloud opacities for each t_p point. Structure is as follows:
            mass_frac[i][j][w][h][k]
            - i: [2 elements] 0 for morning terminator, 1 for evening terminator
            - j: [lat_points elements] latitude point starting from lowest
                 (closeset to lat = -90) to highest (closest to lat = 90)
            - w: [number of wavelength bins] this has to be equivalent to the wavelength
                 structure of prt.
            - h: [number of pressure points] this needs to be equivalent to the pressure
                structure of the gcm.
            - k: 0 for kappa_absorption and 1 for kappa_scatering
        use_bruggemann: bool, optional
            If this flag is set to true, bruggemann mixing is used. This is
            more accurate but takes much longer to calculate.

        Returns
        -------
        wavelengths: xr.DataArray
            Wavelengths of the output spectra in micron
        spectra: xr.DataArray
            transit radius in cm

        """

        # check if data was appropriatly prepared
        if 'transit' not in self.dsi.attrs:
            raise ValueError('To use PrtInterface.calc_transit_spectrum() '
                             'the data needs to be prepared with set_data(..., '
                             'terminator_avg=True, ...).')

        # check if gravity is given, if not use gcm value
        if gravity is None:
            gravity = self.dsi.attrs.get(c["g"]) * 100

        # check if radius is given, if not use gcm value
        if rplanet is None:
            rplanet = self.dsi.attrs.get(c["R_p"]) * 100

        # check if pressure is given, if not use gcm value
        if pressure_0 is None:
            pressure_0 = np.max(self.prt.press)*1e-6

        # check if clouds are wished
        do_clouds = False
        if clouds is not None:
            # check if clouds are possible
            if not is_the_data_cloudy(self.dsi, clouds_only=True):
                raise ValueError('Not all required data to consider clouds within '
                                 'petitRADTRANS are available within the GCM given.')
            # sort out if input was bool or list
            if isinstance(clouds, bool):
                # set cloud flag to input
                do_clouds = clouds
            else:
                do_clouds = True

        # get profile for each latitude and each terminator
        spectra_list = []
        for i, lat in enumerate(self.dsi[c['lat']].values):
            # add morning terminator spectra
            spectra_list.append(self._get_1_transit_spectra([i, lat], [0, -90], mass_frac, gravity,
                                                            mmw, rplanet, pressure_0,
                                                            do_clouds, clouds, use_bruggemann))

            # add evening terminator spectra
            spectra_list.append(self._get_1_transit_spectra([i, lat], [1, 90], mass_frac, gravity,
                                                            mmw, rplanet, pressure_0,
                                                            do_clouds, clouds, use_bruggemann))

        # avarage over all profiles (area avaraged)
        spectra = (np.asarray(spectra_list))**2/len(np.asarray(spectra_list))
        spectra = np.sqrt(np.sum(spectra, axis=0))

        # calcualte wavelengths in micron
        wavelengths = 29979245800.0/self.prt.freq/1e-4

        # return the final spectra
        return wavelengths, spectra

    def _get_1_transit_spectra(self, lat, lon, mass_frac, gravity, mmw, rplanet,
                               pressure_0, do_clouds, clouds, use_bruggemann):
        """
        Calculate the transit spectrum. This function avarages T-p profiles in
        the terminator region and uses poorman equilbriums chemistry.

        Parameters
        ----------
        lat: Tuple
            Tuple containing latitude index and latitude angle.
        lon: Tuple
            Tuple containing longitude index and latitude angle.
        mass_frac: Union[None, List[List[Dict]]]
            The mass fractions for each t_p point. Structure is as follows: mass_frac[i][j]['key']
            - i: [2 elements] 0 for morning terminator, 1 for evening terminator
            - j: [lat_points elements] latitude point starting from lowest
                 (closeset to lat = -90) to highest (closest to lat = 90)
            - dict must contain all mass fractions  elements given as opacity species to prt
        gravity: float
            surface gravity in !cgs!.
        mmw: float or 1D-array
            Mean molecular weight (in atomic units). Will be globally uniform if
            float or horizonatally uniform if 1D.
        rplanet: float
            planet radius in !cm!.
        pressure_0: float
            pressure level of the gravity and radius.
        do_clouds: bool
            Flag if clouds should be included or not
        clouds: Union[np.ndarray[i, j, h, w, k], bool]
            The cloud opacities for each t_p point. Structure is as follows:
            mass_frac[i][j][w][h][k]
            - i: [2 elements] 0 for morning terminator, 1 for evening terminator
            - j: [lat_points elements] latitude point starting from lowest
                 (closeset to lat = -90) to highest (closest to lat = 90)
            - w: [number of wavelength bins] this has to be equivalent to the wavelength
                 structure of prt.
            - h: [number of pressure points] this needs to be equivalent to the pressure
                structure of the gcm.
            - k: 0 for kappa_absorption and 1 for kappa_scatering
        use_bruggemann: bool
            If this flag is set to true, bruggemann mixing is used. This is
            more accurate but takes much longer to calculate.

        Returns
        -------
        prt.transmision object
            petitRADTRAN transmission object

        """
        # get temperature profile
        temp = self.dsi.sel(lat=lat[1], lon=lon[1])[c['T']].values

        # calculate chemistry in mass fractions
        if mass_frac is None:
            chem = self.chemistry.abunds.sel(lat=lat[1], lon=lon[1])
            abus = {}

            # from xarray to dict
            for key in chem.data_vars:
                abus[key] = chem[key].values
            del abus['T']

            # flip the script
            for key in abus:
                abus[key] = abus[key][::-1]
        else:
            abus = mass_frac[0][lat[0]]

        # prepare array for cloud data
        cloud_data = np.zeros((2, len(self.prt.freq), len(self.prt.press)))

        # check if clouds are wished
        if do_clouds:
            # adjust prt for clouds
            self.prt.clouds_mix_opa = types.MethodType(patch_cloud_mix_opa, self.prt)
            self.prt.calc_transm = types.MethodType(patch_calc_transm, self.prt)
            self.prt.delete_clouds = types.MethodType(patch_delete_clouds, self.prt)

            # load or calculate cloud opacities
            if isinstance(clouds, bool):
                # find all volume fractions
                vol_fracs = {}
                for key in self.dsi.keys():
                    if 'ClVf' in key:
                        vol_fracs[key[5:]] = self.dsi.sel(lat=lat[1], lon=lon[1])[key].values[::-1]

                # calculate cloud opacities
                # qabs is the absorption efficiency, qsca the scattering efficiency
                # and csec the cross-section
                qabs, qsca, csec = cloud_opacities(
                    299792458 / self.prt.freq,
                    self.dsi.sel(lat=lat[1], lon=lon[1])['ClDs'].values[::-1]*1e-6,
                    self.dsi.sel(lat=lat[1], lon=lon[1])['ClAb'].values[::-1],
                    self.dsi.sel(lat=lat[1], lon=lon[1])['ClDr'].values[::-1],
                    vol_fracs,
                    use_bruggemann
                    )

                # save opacity kappas
                for k, _ in enumerate(csec):
                    cloud_data[0, :, k] = qabs[k] * csec[k]
                    cloud_data[1, :, k] = qsca[k] * csec[k]
            else:
                # load data from input if given
                cloud_data[0, :, :] = clouds[lon[0], lat[0], :, :, 0]
                cloud_data[1, :, :] = clouds[lon[0], lat[0], :, :, 1]

        # check if all opacity species are in mass fractions
        for key in self.prt.line_species:
            if key not in abus:
                # if a species is missing, check for partial matches
                for k_abus in abus:
                    if k_abus in key:
                        abus[key] = abus[k_abus]
                        break
                else:
                    # if a species can not be found fully or partally, give error
                    raise ValueError('Opacity species ' + key + ' missing in mass fraction input')

        # calculate transmission spectra for the morning terminator of current lat point
        if do_clouds:
            self.prt.calc_transm(temp[::-1], abus, gravity, mmw*np.ones_like(temp),
                                 R_pl=rplanet, P0_bar=pressure_0, clouds=cloud_data)
        else:
            self.prt.calc_transm(temp[::-1], abus, gravity, mmw*np.ones_like(temp),
                                R_pl=rplanet, P0_bar=pressure_0)

        # add the spectra to the list
        return self.prt.transm_rad

    def _get_stellar_spec(self, wlen, t_star):
        """
        Helperfunction from petitRADTRANS that calculates the stellar spectrum
        from the phoenix spectrum
        """
        from petitRADTRANS.nat_cst import get_PHOENIX_spec_rad

        if t_star is not None:
            spec, _ = get_PHOENIX_spec_rad(t_star)
            stellar_intensity = np.interp(wlen * 1e-4, spec[:, 0], spec[:, 1])
            return stellar_intensity

        raise ValueError('Tstar is need to define a stellar spectra.')


class PACInterface(Interface):
    """
    Interface GCM data to a 1D or 2D PAC chemistry simulation.
    (Note: PAC is not needed for these routines to work.)
    """

    def __init__(self, tools, pac_dim, dsi=None):
        """
        Constructs the Interface and links to either 1D or 2D PAC.

        Parameters
        ----------
        tools: GCMTools Object
            A GCMTools object that is linked to the interface
        pac_dim: int
            Type of PAC simulation '1' or '2' D
        dsi: DataSet, optional
            Shortcut to set the dataset to the one that is required.
            (Otherwise, use the set_data method.)
        """
        super().__init__(tools)

        if not (pac_dim == 1 or pac_dim == 2):
            raise ValueError("Please enter a valid PAC dimension: '1' or '2' \
                             for 1D or pseudo-2D.")
        self.dim = pac_dim

        if dsi is not None:
            self.dsi = dsi

    def write_inputfile(self, destination, kwargs_1D={}, kwargs_2D={}):
        """
        Transform the given GCM dataset to a pseudo-2D PAC
        input file.

        Parameters
        ----------
        destination: str
            Path of the directory where the input file should be written to.
        kwargs_1D: dict, optional
            A dictionary containing all of the necessary keyword arguments for
            a 1D PAC input file.
        kwargs_2D: dict, optional
            A dictionary containing all of the necessary keyword arguments for
            a pseudo-2D PAC input file.
        """
        import os

        # check if data is present
        if self.dsi is None:
            raise RuntimeError("First select the required dataset with the \
                                set_data method.")
        # check if the path exists
        if not os.path.isdir(destination):
            raise OSError("The given destination directory does not exist:\n" +
                          destination)

        # call the right method
        if self.dim == 1:
            self._write_1Dpac_inpfile(self.dsi, destination, **kwargs_1D)
        elif self.dim == 2:
            self._write_2Dpac_inpfile(self.dsi, destination, **kwargs_2D)

    def _write_1Dpac_inpfile(self, dsi, destination,
                            spec_file='', zab_file='', reac_file='', therm_file='',
                            eddy_file='', star_file='',
                            pressure_bot=None, pressure_top=None, np=None,
                            R_star=1.0, R_planet=11.2, M_planet=317.8, a=0.01,
                            zenith_angle=48.0, albedo=0.0, ipho_file='',
                            use_mixing=True, use_photo=True, model_name=None):
        """
        Write a 1D PAC-input file based on the given GCM dataset.

        Parameters
        ----------
        dsi: DataSet
            GCM dataset that forms the basis of the input file.
        destination: str
            Path of the directory where the input file should be written to.
        spec_file: string, optional
            The name of the file containing abundances and species.
            If no name is given, a placeholder name is written.
        zab_file: string, optional
            The name of the zab file (generated with ACE code).
            If no name is given, a placeholder name is written.
        reac_file: string, optional
            The name of the reac file, containing all chemical reactions.
            If no name is given, a placeholder name is written.
        therm_file: string, optional
            The name of the file containing the NASA thermodynamic data.
            If no name is given, a placeholder name is written.
        eddy_file: string, optional
            The name of the file containg eddy diffusion coefficients.
            If no name is given, a placeholder name is written.
        star_file: string, optional
            The name of the file containing the stellar spectrum.
            If no name is given, a placeholder name is written.
        pressure_bot: float, optional
            Bottom/maximum pressure of the grid (in bar).
            If no number is given, the bottom pressure of the dataset is used.
        pressure_top: float, optional
            Top/minimum pressure of the grid (in bar).
            if no number is given, the top pressure of the dataset is used.
        np: int, optional
            Number of pressures in the vertical grid.
            If no number is given, it is the same as that of the dataset.
        R_star: float, optional
            Stellar radius in solar radii. Default: 1 solar radius.
        R_planet: float, optional
            Planetary radius in Earth radii. Default: 11.2 (radius of Jupiter).
        M_planet: float, optional
            Planetary mass in Earth masses. Default: 317.8 (mass of Jupiter).
        a: float, optional
            Semi-major axis in AU. Default: 0.01.
        zenith_angle: float, optional
            Zenith angle for the irradiation in degrees. Default: 48 deg.
        albedo: float, optional
            Surface albedo of the planet. Default: 0.0 (no reflection).
        ipho_file: string, optional
            File containing all of the individual photochemical dissociations.
            If no name is given, a placeholder name is written.
        use_mixing: boolean, optional
            Whether to use vertical mixing or not. Default: True.
        use_photo: boolean, optional
            Whether to use photochemistry/dissociations or not. Default: True
        model_name: string, optional
            The name of the chemistry model: *model_name*.inp.
            If no name is given, the 'tag' of the dataset is used.
        """
        # Checks and assigning default values
        if model_name is None:
            model_name = dsi.tag

        if not spec_file:
            spec_file = 'template.spec'
        if not zab_file:
            zab_file = model_name + '.zab'
        if not reac_file:
            reac_file = 'template.reac'
        if not therm_file:
            therm_file = 'template.therm'
        if not eddy_file:
            eddy_file = model_name + '.eddy'
        if not star_file:
            star_file = 'template.star'
        if not ipho_file:
            ipho_file = 'template.ipho'

        if pressure_bot is None:
            pressure_bot = max(dsi.Z.values)
        if pressure_top is None:
            pressure_top = min(dsi.Z.values)
        if dsi.p_unit == 'Pa':  # convert to bar if needed
            pressure_bot = pressure_bot / 1e5
            pressure_top = pressure_top / 1e5
        if np is None:
            np = len(dsi.Z.values)

        # Write all parameters in the required format as output file
        with open(destination +'/'+model_name+'.inp', 'w') as f:
            f.write(spec_file + (35-len(spec_file))*' ' + '! species file\n')
            f.write(zab_file + (35-len(zab_file))*' ' + '! z,p,T (initial abundances) file\n')
            f.write('2  0                               ! No reaction files, (0/1) not/use iondip treatment\n')
            f.write(reac_file + (35-len(reac_file))*' ' + '! reaction file\n')
            f.write(therm_file + (35-len(therm_file))*' ' + '! reaction file\n')
            f.write(eddy_file + ' 0' + (33-len(eddy_file))*' ' + '! eddy diffusion coefficient profile file\n')
            f.write(star_file + (35-len(star_file))*' ' + '! stellar spectrum file\n')
            f.write(ipho_file + (35-len(ipho_file))*' ' + '! photo cross sections info file\n')
            f.write(f'{pressure_bot:.1f}d0  {pressure_top:.4e}  {np}            ! pressure [bar] at bottom/top, No heights\n')
            f.write(f'{R_star:.3f}d0                            ! star radius [R(Sun)]\n')
            f.write(f'{R_planet:.3f}d0  {M_planet:.2f}d0                 ! planet radius [R(Earth)], planet mass [m(Earth)]\n')
            f.write(f'{a:.5f}d0  {zenith_angle:.1f}d0  {albedo:.1f}d0           ! orbital distance [AU], zenith angle [deg], surface albedo\n')
            f.write(f'{int(use_mixing)}  {int(use_photo)}                               ! deactivate/activate (0/1) diffusion and photochemistry\n')
            f.write('2                                  ! numerical method (1/2/3) (-1 is solid body rotation)\n')

        wrt.write_status('INFO', 'File written: '+destination+'/'+model_name+'.inp')

    def _write_2Dpac_inpfile(self, dsi, destination,
                            spec_file='', zab_file='', reac_file='', therm_file='',
                            eddy_file='', star_file='', lpt_file='',
                            pressure_bot=None, pressure_top=None, np=None,
                            R_star=1.0, R_planet=11.2, M_planet=317.8, a=0.01,
                            zenith_angle=48.0, albedo=0.0,
                            nlon=90, nrot=30, ipho_file='', use_2D_eddy=False,
                            use_mixing=True, use_photo=True, model_name=None):
        """
        Write a pseudo-2D PAC-input file based on the given GCM dataset.

        Parameters
        ----------
        dsi: DataSet
            GCM dataset that forms the basis of the input file.
        destination: str
            Path of the directory where the input file should be written to.
        spec_file: string, optional
            The name of the file containing abundances and species.
            If no name is given, a placeholder name is written.
        zab_file: string, optional
            The name of the zab file (generated with ACE code).
            If no name is given, a placeholder name is written.
        reac_file: string, optional
            The name of the reac file, containing all chemical reactions.
            If no name is given, a placeholder name is written.
        therm_file: string, optional
            The name of the file containing the NASA thermodynamic data.
            If no name is given, a placeholder name is written.
        eddy_file: string, optional
            The name of the file containg eddy diffusion coefficients.
            If no name is given, a placeholder name is written.
        star_file: string, optional
            The name of the file containing the stellar spectrum.
            If no name is given, a placeholder name is written.
        lpt_file: string, optional
            The name of the lpt-file, containing longitude-pressure-temperature.
            If no name is given, a placeholder name is written.
        pressure_bot: float, optional
            Bottom/maximum pressure of the grid (in bar).
            If no number is given, the bottom pressure of the dataset is used.
        pressure_top: float, optional
            Top/minimum pressure of the grid (in bar).
            if no number is given, the top pressure of the dataset is used.
        np: int, optional
            Number of pressures in the vertical grid.
            If no number is given, it is the same as that of the dataset.
        R_star: float, optional
            Stellar radius in solar radii. Default: 1 solar radius.
        R_planet: float, optional
            Planetary radius in Earth radii. Default: 11.2 (radius of Jupiter).
        M_planet: float, optional
            Planetary mass in Earth masses. Default: 317.8 (mass of Jupiter).
        a: float, optional
            Semi-major axis in AU. Default: 0.01.
        zenith_angle: float, optional
            Zenith angle for the irradiation in degrees. Default: 48 deg.
        albedo: float, optional
            Surface albedo of the planet. Default: 0.0 (no reflection).
        nlon: int, optional
            Number of longitudinal cells in the chemistry code. Default: 90.
        nrot: int, optional
            Number of rotations the code should integrate for. Default: 30.
        ipho_file: string, optional
            File containing all of the individual photochemical dissociations.
            If no name is given, a placeholder name is written.
        use2Deddy: boolean, optional
            A 1D (False; default) or 2D (True) eddy diffusion profile is used.
        use_mixing: boolean, optional
            Whether to use vertical mixing or not. Default: True.
        use_photo: boolean, optional
            Whether to use photochemistry/dissociations or not. Default: True
        model_name: string, optional
            The name of the chemistry model: *model_name*.inp.
            If no name is given, the 'tag' of the dataset is used.
        """
        # Checks and assigning default values
        if model_name is None:
            model_name = dsi.tag

        if not spec_file:
            spec_file = 'template.spec'
        if not zab_file:
            zab_file = model_name + '.zab'
        if not reac_file:
            reac_file = 'template.reac'
        if not therm_file:
            therm_file = 'template.therm'
        if not eddy_file:
            eddy_file = model_name + '.eddy'
        if not star_file:
            star_file = 'template.star'
        if not lpt_file:
            lpt_file = model_name + '.lpt'
        if not ipho_file:
            ipho_file = 'template.ipho'

        if pressure_bot is None:
            pressure_bot = max(dsi.Z.values)
        if pressure_top is None:
            pressure_top = min(dsi.Z.values)
        if dsi.p_unit == 'Pa':  # convert to bar if needed
            pressure_bot = pressure_bot / 1e5
            pressure_top = pressure_top / 1e5
        if np is None:
            np = len(dsi.Z.values)

        # Write all parameters in the required format as output file
        with open(destination+'/'+model_name+'.inp', 'w') as f:
            f.write(spec_file + (35-len(spec_file))*' ' + '! species file\n')
            f.write(zab_file + (35-len(zab_file))*' ' + '! z,p,T (initial abundances) file\n')
            f.write('2  0                               ! No reaction files, (0/1) not/use iondip treatment\n')
            f.write(reac_file + (35-len(reac_file))*' ' + '! reaction file\n')
            f.write(therm_file + (35-len(therm_file))*' ' + '! reaction file\n')
            f.write(eddy_file + f' {int(use_2D_eddy)}' + (33-len(eddy_file))*' ' + '! eddy diffusion coefficient profile file\n')
            f.write(star_file + (35-len(star_file))*' ' + '! stellar spectrum file\n')
            f.write(ipho_file + (35-len(spec_file))*' ' + '! photo cross sections info file\n')
            f.write(f'{pressure_bot:.1f}d0  {pressure_top:.4e}  {np}            ! pressure [bar] at bottom/top, No heights\n')
            f.write(f'{R_star:.3f}d0                            ! star radius [R(Sun)]\n')
            f.write(f'{R_planet:.3f}d0  {M_planet:.2f}d0                 ! planet radius [R(Earth)], planet mass [m(Earth)]\n')
            f.write(f'{a:.5f}d0  {zenith_angle:.1f}d0  {albedo:.1f}d0           ! orbital distance [AU], zenith angle [deg], surface albedo\n')
            f.write(f'{int(use_mixing)}  {int(use_photo)}                               ! deactivate/activate (0/1) diffusion and photochemistry\n')
            f.write('-1                                 ! numerical method (1/2/3) (-1 is solid body rotation)\n')
            f.write(lpt_file + (35-len(lpt_file))*' ' + '! sbr: longitude-pressure-temperature structure file\n')
            f.write(f'{nlon} {nrot}                              ! sbr: No longitudes per period, No rotation periods')

        wrt.write_status('INFO', 'File written: '+destination+'/'+model_name+'.inp')


    def generate_lptfile(self, destination, jet_speed=None, set_min_temp=None,
                         eps=20, tave=False, kwargs_thermosphere={},
                         plot_input=False, model_name=None):
        """
        Generate from the selected dataset a pseudo-2D PAC (Agundez+2014)
        'lpt' input file (longitude [in units of pi], pressure [in
        bar], temperature [in K]).

        The zonal wind speed for the solid-body rotation of the chemical
        atmosphere is extracted from the dataset. It is also possible to specify a value.

        Parameters
        ----------
        destination: str
            Path of the directory where the input file should be written to.
        jet_speed: float, optional
            Zonal wind value (in m/s) that is used in the solid-body rotation
            mode of PAC. If no value is given, this parameter is extracted from
            the GCM dataset.
        set_min_temp: float, optional
            If given, the atmosphere temperature is set to this value whenever
            it would actually be lower.
        eps: float, optional
            Epsilon value (in degrees latitude). The meridional mean is
            calculated between -eps and +20 around the equator.
        tave: boolean, optional
            Flag whether the regular snapshot data should be used (default), or
            the time averaged data.
        kwargs_thermosphere: dict, optional
            Optional set of arguments used to extend the atmosphere upwards.
        plot_input: boolean, optional
            Set to true if a plot of the input temperature data should be made,
            in the same directory as the lpt-file.
        model_name: string, optional
            The name of the chemistry model: *model_name*.lpt.
            If no name is given, the 'tag' of the dataset is used.
        """
        # Check if data is present
        if self.dsi is None:
            raise RuntimeError("First select the required dataset with the \
                                set_data method.")
        # Check if the path exists
        import os
        if not os.path.isdir(destination):
            raise OSError("The given destination directory does not exist:\n" +
                          destination)
        # Checks and assign model name
        if model_name is None:
            model_name = self.dsi.tag
        # Assign temperature key
        if tave:
            t_key = 'Ttave'
        else:
            t_key = 'T'

        # If needed, derive the zonal wind speed from the dataset
        if not jet_speed:
            jet_speed = self._extract_jet_speed(self.dsi, eps=eps, tave=tave)

        # If required, extend temperatures upward with a thermosphere
        if kwargs_thermosphere:
            from ..utils import manipulations as man
            T = man.m_extend_upward(self.dsi[t_key], **kwargs_thermosphere)
        # If not, just work with the native temperature DataArray
        else:
            T = self.dsi[t_key]

        # Extract pressures
        if self.dsi.p_unit == 'Pa':  # convert to bar if needed
            p_bar = [ip/1e5 for ip in T.Z.values]
        elif self.dsi.p_unit == 'bar':
            p_bar = T.Z.values

        # Define longitude sampling
        lon_grid = np.linspace(-180, 178, 180)
        # Again, but units of pi and start at substellar point: 0 (pi) = 2 (pi)
        lon_grid_out = np.linspace(0, 2, 181)

        # Interpolate temperature grid to new lon coordinates
        # NOTE: Careful with extrapolation at the edges of the grid!
        #       Cyclic interpolation would probably be more robust...
        T = T.interp(lon=lon_grid, kwargs={'fill_value':'extrapolate'})
        # Meridional mean of the equatorial region
        T = T.sel(lat=slice(-eps,eps))
        weights = np.cos(np.deg2rad(T.lat))
        T = T.weighted(weights).mean(dim='lat')
        # 'Roll' the data so that substellar is left and antistellar is center
        T = T.roll(lon=-90, roll_coords=True)
        T = T.values

        # Restrict if a minimum temperature is given
        if set_min_temp:
            T[T < set_min_temp] = set_min_temp

        # Writing the lpt-file
        with open(destination+'/'+model_name+'.lpt', 'w') as f:
            # header, zonal wind and grid info
            f.write(f'  {jet_speed/1000:1.3f}     ! velocity [km/s]\n')
            f.write(f'  {len(lon_grid_out)}       ! number of longitudes\n')
            f.write(f'  {len(p_bar)}        ! number of pressures\n')
            # pressures
            f.write('! next line lists pressures [bar]\n')
            for ip in p_bar:
                f.write(f'     {ip:1.5E}')
            f.write('\n')
            # temperatures per longitude
            f.write('! longitude[pi]          Tk[K]          ...\n')
            for i in range(0, len(lon_grid)):
                f.write(f'     {lon_grid_out[i]:1.4f}     ')
                for k in range(0, len(p_bar)):
                    f.write(f'      {T[k,i]:4.2f}')
                f.write('\n')
            # repeat first entry (substellar point)
            f.write(f'     {lon_grid_out[-1]:1.4f}     ')
            for k in range(0, len(p_bar)):
                f.write(f'      {T[k,0]:4.2f}')

        wrt.write_status('INFO', 'File written: '+destination+'/'+model_name+'.lpt')

        # If required, also save a plot of the input data
        if plot_input:
            import matplotlib.pyplot as plt
            fig = plt.figure('Lpt plot')
            ax = plt.gca()
            image = ax.pcolormesh(lon_grid_out, p_bar, T, cmap='inferno',
                        linewidth=0, rasterized=True) # these options reduce
                                                      # visual grid glitches
            image.set_edgecolor('face')
            cb = plt.colorbar(image)
            # in-plot wind info
            ax.text(lon_grid_out[2], p_bar[5], f'U = {jet_speed/1000:.1f} km/s',
                    fontsize=12)
            ax.set_yscale('log')
            ax.invert_yaxis()
            ax.set_xlabel(r'Longitude from substellar point ($\pi$)', fontsize=14)
            ax.set_ylabel('Pressure (bar)', fontsize=14)
            cb.set_label('Temperature (K)', fontsize=14)
            plt.savefig(destination+f'/lpt_plot.pdf', format='pdf')
            plt.close(fig)


    def _extract_jet_speed(self, dsi, eps=20, tave=False):
        """Function to extract the jet stream speed in m/s as a single value
        from the given dataset.

        Parameters
        ----------
        dsi: DataSet
            The GCM from which the wind speed is to be extracted.
        eps: float, optional
            Epsilon value (in degrees latitude). The meridional mean is
            calculated between -eps and +20 around the equator.
        tave: boolean, optional
            Flag whether the regular snapshot data should be used (default), or
            the time averaged data.
        Returns
        -------
        jet_speed: float
            A single value for the zonal mean equatorial jet speed (in m/s).
        """
        # Assign the correct variable key
        if tave:
            u_key = 'uVeltave'
        else:
            u_key = 'U'

        weights = np.cos(np.deg2rad(dsi.lat))   # weight with cos(latitude)
        U = dsi[u_key].sel(lat=slice(-eps,eps)) # take latitudinal band from -20 to +20 degrees (around equator)
        U = U.weighted(weights).mean(dim='lat') # average latitudinal band meridonally

        # average over pressure levels (10 bar and up)
        if dsi.p_unit == 'Pa':  # convert to bar if needed
            p_max = 1e6
        elif dsi.p_unit == 'bar':
            p_max = 10
        U = U.sel(Z=slice(p_max, None)).mean(dim='Z')
        jet_speed = U.mean(dim='lon')           # average end result zonally
        jet_speed = jet_speed.values

        # if result turns out to be negative, raise a warning
        if jet_speed < 0:
            msg = f'WARNING: zonal wind speed {zonalwindspeed} is negative!\n'
            msg += 'Setting wind speed to 10 m/s.'
            wrt.write_message(msg, color='WARN')
            zonalwindspeed = 10

        return jet_speed

    def generate_aptfiles(self, destination, lons=[], eps=20, tave=False,
                          kwargs_thermosphere={}, plot_input=False,
                          model_name=None):
        """Function to write apt-files (altitude [in km], pressure [in bar],
        temperature [in K]), to use as input in the ACE chemical equilibrium
        code.
        A text-file summarizing all longitudes is also written.

        Parameters
        ----------
        destination: str
            Path of the directory where the input file should be written to.
        lons: [float], optional
            List of the longitudes of the vertical columns of this dataset that
            need to be sampled and written to apt-files.
            If multiple are given, multiple apt-files will be written.
            If none are given, the substellar point is selected.
        eps: float, optional
            Epsilon value (in degrees latitude). The meridional mean is
            calculated between -eps and +20 around the equator.
        tave: boolean, optional
            Flag whether the regular snapshot data should be used (default), or
            the time averaged data.
        kwargs_thermosphere: dict, optional
            Optional set of arguments used to extend the atmosphere upwards.
        plot_input: boolean, optional
            Set to true if a plot of the input temperature data should be made,
            in the same directory as the apt-file.
        model_name: string, optional
            The name of the chemistry model: *model_name*_*lon*.apt.
            If no name is given, the 'tag' of the dataset is used.
        """
        # Check if data is present
        if self.dsi is None:
            raise RuntimeError("First select the required dataset with the \
                                set_data method.")
        # Check if the path exists
        import os
        if not os.path.isdir(destination):
            raise OSError("The given destination directory does not exist:\n" +
                          destination)
        # Checks and assign model name
        if model_name is None:
            model_name = self.dsi.tag
        # Assign temperature key
        if tave:
            t_key = 'Ttave'
        else:
            t_key = 'T'
        # If a previous list of apt-files is still present, remove it.
        try:
            os.remove(destination + '/apt_list.txt')
        except OSError:
            pass
        # If no longitude is given, just take the substellar point at 0 degrees
        if not list(lons):
            lons=[0]
        # If an overview plot of the input is needed, prepare it here.
        if plot_input:
            import matplotlib
            import matplotlib.pyplot as plt
            fig_many = plt.figure()
            ax_many = plt.gca()
            norm = matplotlib.colors.Normalize(vmin=0, vmax=180)
            cmap = matplotlib.cm.get_cmap('cividis')

        # If required, extend temperatures upward with a thermosphere
        if kwargs_thermosphere:
            from ..utils import manipulations as man
            Tsource = man.m_extend_upward(self.dsi[t_key], **kwargs_thermosphere)
        else:
            Tsource = self.dsi[t_key]

        # For each longitude given...
        for lon in lons:
            # Extract the corresponding longitude
            T = Tsource.sel(lon=lon, method='nearest').sel(lat=slice(-eps, eps))
            w = np.cos(np.deg2rad(self.dsi.lat))     # weighted with cos(lat)
            T = T.weighted(w).mean(dim='lat')        # meridional average

            # Calculate height based on hydrostatic equilibrium, and renormalize
            # so that 1 bar <--> 0 m
            R_p = self.dsi.R_p            # radius in m
            G = 6.674e-11                 # gravitational const in m3 kg−1 s−2
            M_p = self.dsi.g * R_p**2 / G # mass in kg
            R_spec = self.dsi.R           # specific gas constant in J kg-1 K-1
            p_ref = 1e5                   # reference pressure in pascal (1 bar)
            # pressures should be in bar, except for the altitude-routine
            if self.dsi.p_unit == 'bar':
                p = T.Z.values
                p_pascal = [ip*1e5 for ip in p]
            elif self.dsi.p_unit == 'Pa':
                p_pascal = T.Z.values
                p = [ip/1e5 for ip in p_pascal]

            alt = self._p_to_alt(p_pascal, T, M_p, R_p, R_spec, p_ref)

            # Construct full filename:
            full_path = destination + '/' + model_name + f'_{int(lon)}.apt'

            # Write the file:
            # should be written in order of decreasing pressures:
            if p[1] < p[0]:
                index_order = range(0, len(p))             # write in same order
            elif p[1] > p[0]:
                index_order = range(len(p)-1, -1, -1)   # write in reverse order
            with open(full_path, 'w') as f:
                f.write('! altitude[km]   pressure[bar]   temperature[K]\n')
                for k in index_order:
                    line = '  ' + '{:4.4f}'.format(alt[k]/1000) + '   ' + \
                            '  ' + '{:1.4E}'.format(p[k]) + '  ' + \
                            '  ' + '{:4.2f}'.format(T.values[k]) + '\n'
                    f.write(line)
            wrt.write_status('INFO', 'File written: ' + full_path)

            # Remember the apt-file in a list for easy submission
            with open(destination+'apt_list.txt', 'a') as listfile:
                listfile.write(model_name + f'_{int(lon)}\n')

            # If required, a control plot is saved to show the input data
            if plot_input:
                mycol = cmap(norm(180-abs(lon)))
                if lon > 0:
                    ls = '--'
                else:
                    ls = '-'
                ax_many.plot(T.values, p, ls, color=mycol, linewidth=2)

        if plot_input:
            ax_many.set_yscale('log')
            ax_many.invert_yaxis()
            ax_many.set_xlabel('Temperature (K)', fontsize=16)
            ax_many.set_ylabel('Pressure (bar)', fontsize=16)
            plt.savefig(destination+f'PT_profile_all.pdf', format='pdf')
            plt.close(fig_many)


    def _p_to_alt(self, p, T, M_p, R_p, R_spec=3589, p_ref=1.e7):
        """ Compute the corresponding height (in m) for a given list of
        pressures, assuming a hydrostatic atmosphere. The gravity changes as a
        function of the vertical layers (R + dz), but is assumed constant within
        one layer. The surface gravity (g = G M/ R^2) is assumed to be at the
        reference pressure. The algorithm works in 4 parts:
            1) a pivot near the reference pressure (= surface gravity) is found
            2) integrated from pivot upward
            3) integrated from pivot downward
            4) small correction because the pivot is not exactly at p_ref

        Parameters
        ----------
        p: [float]
            Pressures (in pascal!).
        T: [float]
            Corresponding temperatures in kelvin.
        M_p: float
            Planet mass in kg.
        R_p: float
            Planet radius in m.
        R_spec: float, optional
            Specific gas constant of the atmosphere (R_spec=2/7*cp for diatomic gas).
        p_ref: float, optional
            Reference pressure in pascal.

        Returns
        -------
        alt_scaled: [float]
            List of altitudes (in meter) corresponding to the given pressures.
        """
        from scipy import interpolate

        # Find index of the reference pressure (gravity = 0)
        for i in range(0, len(p)):
            if p[i] < p_ref:
                pivot = i
                break

        # Initialize height profile
        z = np.zeros(len(p))

        # Reconstructing the height with numerical integration (upper part)
        for i in range(pivot, len(p)):
            if i == pivot:  # first iteration
                gl = 6.674e-11 * M_p / (R_p)**2     # surface gravity

                Tp = interpolate.interp1d(p, T, 'quadratic')    # find T(p_ref)
                T_pref = Tp(p_ref)
                Tl = 0.5*(T[i] + T_pref)    # average temperature in this layer
                Hl = R_spec * Tl / gl       # average scale height in this layer

                # calculate first height as layer thickness between p(pivot) and
                # p_ref
                z[i] = Hl*math.log(p_ref/p[i])

            else: # calculate height for i
                gl = 6.674e-11 * M_p / (R_p + z[i-1])**2    # gravity at bottom
                                                            # edge of this layer

                Tl = 0.5*(T[i] + T[i-1])    # average temperature in this layer
                Hl = R_spec * Tl / gl       # average scale height in this layer

                # calculate new height as previous + layer thickness
                z[i] = z[i-1] + Hl * math.log(p[i-1]/p[i])

        # Reconstructing the height with numerical integration (bottom part)
        for i in range(pivot-1, -1, -1):
            gl = 6.674e-11 * M_p / (R_p + z[i+1])**2    # gravity at top
                                                        # edge of this layer

            Tl = 0.5*(T[i] + T[i+1])    # average temperature in this layer
            Hl = R_spec * Tl / gl       # average scale height in this layer

            # calculate new height as above - new layer thickness
            z[i] = z[i+1] - Hl * math.log(p[i]/p[i+1])

        # Apply (small) offset to the reference pressure
        hf = interpolate.interp1d(p, z, 'quadratic')
        offset = hf(p_ref)
        alt_scaled = [thisz - offset for thisz in z]

        return alt_scaled
