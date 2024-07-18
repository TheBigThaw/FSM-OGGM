"""
Mass balance class FSM - OGGM coupling
"""
# External libraries
import logging
import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import warnings
import salem
from oggm.core.massbalance import MassBalanceModel
from oggm.utils import ncDataset
from oggm import cfg, utils
from oggm.exceptions import (InvalidParamsError)
import xarray as xr
import glob
from functools import partial
import FSM

cfg.add_to_basenames('WFDE5_Hintereisferner_1980-2019',
                     'WFDE5_Hintereisferner_1980-2019.nc',
                     'FSM ready climate file')
cfg.add_to_basenames('climate_historical_fsm',
                     'climate_historical_fsm.nc',
                     'FSM ready climate file')

# Module logger
log = logging.getLogger(__name__)


def _preprocess(x, i, j):
    """
    Preprocess to pass to xarray.open_mfdataset()
    :param x: Xarray Dataset
    :param i: coordinate for longitude
    :param j: coordinate for latitude
    :return: a selected xarray Dataset for two specific coordinates
    """
    return x.isel(lat=j, lon=i)


@utils.entity_task(log, writes=['climate_historical_fsm'])
def process_wfde5_data(gdir,
                       y0=None,
                       y1=None):
    """
    Process WFDE5 meteorological variables and put them in an FSM ready format
    per glacier directory
    """

    output_file_name = gdir.get_filepath('climate_historical_fsm')
    if os.path.exists(output_file_name):
        os.remove(output_file_name)

    lat = gdir.cenlat
    lon = gdir.cenlon

    if y0 is None:
        y0 = '1979'
    if y1 is None:
        y1 = '2019'

    # nearest point on global 0.5 degree grid
    i = round(2 * (179.75 + lon))
    j = round(2 * (89.750 + lat))

    # location and height of reference pixel
    fpath = os.path.join(cfg.PATHS['climate_file'], 'ASurf_WFDE5_CRU_v2.0.nc')
    df = xr.open_dataset(fpath)
    ref_pix_lat = df['lat'][j].values
    ref_pix_lon = df['lon'][i].values
    ref_hgt = df['ASurf'][j, i].values

    paths_files = sorted(glob.glob(os.path.join(cfg.PATHS['climate_file'],
                                                '**/*' + '_WFDE5_CRU_' + '*_v2.0.nc')))

    partial_func = partial(_preprocess, i=i, j=j)

    ds = xr.open_mfdataset(paths_files,
                           combine='nested',
                           concat_dim="time",
                           preprocess=partial_func)

    ds.attrs = {'author': 'Beatriz Recinos and Richard Essery',
                'author_info': 'Big thaw OGGM-FSM',
                'ref_hgt': ref_hgt,
                'ref_pix_lat': ref_pix_lat,
                'ref_pix_lon': ref_pix_lon,
                'climate_source': 'WFDE5',
                'yr_0': y0,
                'yr_1': y1}

    ds.load().to_netcdf(gdir.get_filepath('climate_historical_fsm'))


class FactorialSnowpackModel(MassBalanceModel):
    def __init__(self,
                 gdir,
                 filename='climate_historical_fsm',
                 input_filesuffix='',
                 mb=0.):
        super(FactorialSnowpackModel, self).__init__()
        self.hemisphere = 'nh'
        self.valid_bounds = [-2e4, 2e4]  # in m

        # FSM layers
        self.zmin = 2507  # Centre of lowest elevation band (m)
        self.zmax = 3739  # Centre of highest elevation band (m)
        self.Nbnd = 10  # Number of elevation bands
        self.zbnd = self.zmin + (np.arange(self.Nbnd) + 0.5) * (
                    self.zmax - self.zmin) / self.Nbnd  # Elevations of bands (m)
        self.Dmin = np.array([0.1, 0.2, 0.4], 'f')  # Minimum snow layer thicknesses (m)
        self.Nsmx = len(self.Dmin)  # Maximum number of snow layers
        self.Dice = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0], 'f')  # Ice layer thicknesses (m)
        self.Nice = len(self.Dice)

        # FSM state variables
        self.Tm = 273.15
        self.albs = np.full(self.Nbnd, 0.8, 'f')  # Snow albedo
        self.Nsnw = np.zeros(self.Nbnd, 'i')  # Number of snow layers
        self.Tsrf = np.full(self.Nbnd, self.Tm, 'f')  # Surface temperature (K)
        self.Dsnw = np.zeros((self.Nsmx, self.Nbnd), 'f', order='F')  # Snow layer thicknesses (m)
        self.Sice = np.zeros((self.Nsmx, self.Nbnd), 'f', order='F')  # Ice content of snow layers (kg/m^2)
        self.Sliq = np.zeros((self.Nsmx, self.Nbnd), 'f', order='F')  # Liquid content of snow layers (kg/m^2)
        self.Tice = np.full((self.Nice, self.Nbnd), self.Tm, 'f', order='F')  # Ice layer temperatures (K)
        self.Tsnw = np.full((self.Nsmx, self.Nbnd), self.Tm, 'f', order='F')  # Snow layer temperatures (K)

        # Read climate file
        fpath = gdir.get_filepath(filename, filesuffix=input_filesuffix)
        if not os.path.exists(fpath):
            raise FileNotFoundError('run process_wfde5_data')

        with ncDataset(fpath, mode='r') as nc:
            # time
            self.LW = nc.variables['LWdown'][:]
            self.Ps = nc.variables['PSurf'][:]
            self.Qa = nc.variables['Qair'][:]
            self.Rf = nc.variables['Rainf'][:]
            self.Sf = nc.variables['Snowf'][:]
            self.SW = nc.variables['SWdown'][:]
            self.Ta = nc.variables['Tair'][:]
            self.Ua = nc.variables['Wind'][:]
            self.time = nc.variables['time'][:]
            self.Ntim = int(len(self.time))
            self.zref = nc.getncattr('ref_hgt')
            self.dz = self.zbnd - self.zref
        self._mb = mb

    def get_annual_mb(self):
        mb = FSM.fsmpy(self.Dice, self.Dmin, self.dz, self.LW, self.Ps,
                       self.Qa, self.Rf, self.Sf, self.SW, self.Ta,
                       self.Ua, self.albs, self.Dsnw, self.Nsnw, self.Sice,
                       self.Sliq, self.Tice, self.Tsnw, self.Tsrf, self.Nbnd,
                       self.Nice, self.Nsmx, self.Ntim)
        return mb
