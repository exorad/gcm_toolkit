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
        self._set_data_common(time, tag=tag, regrid_lowres=regrid_lowres)

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
        dsi = self.tools.get_one_model(tag).sel(time=time)

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
        self._set_data_common(time, tag=tag, regrid_lowres=regrid_lowres)

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

        raise ValueError("Tstar is need to define a stellar spectra.")
