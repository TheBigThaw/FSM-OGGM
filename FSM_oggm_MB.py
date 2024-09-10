"""
Mass balance class FSM - OGGM coupling
"""
# External libraries
import logging
import multiprocessing
import os
import numpy as np
from datetime import datetime
from oggm.core.massbalance import MassBalanceModel
from oggm.utils import ncDataset
from oggm.cfg import SEC_IN_YEAR
import netCDF4
from oggm import cfg, utils
from oggm import entity_task
from scipy.interpolate import interp1d
import xarray as xr
import glob
import re
from IPython import embed
from functools import partial
from progressbar import ProgressBar, Percentage, Bar
import FSM

cfg.add_to_basenames('WFDE5_Hintereisferner_1980-2019',
                     'WFDE5_Hintereisferner_1980-2019.nc',
                     'FSM ready climate file')
cfg.add_to_basenames('climate_historical_fsm',
                     'climate_historical_fsm.nc',
                     'FSM ready climate file')

# Module logger
log = logging.getLogger(__name__)


def find_files_per_var(var='', y0=None, y1=None):
    """
    Find file paths per variable in climate dir
    :param var: Climate variable to find
    :param y0: (optional) if None picks 1980
    :param y1: (optional) if None picks 2019
    :return: paths with file names matching the variable and the years selected

    """

    if y0 is None:
        y0 = '1980'
    if y1 is None:
        y1 = '2019'

    if var in ('Rainf', 'Snowf'):
        paths_files = sorted(glob.glob(os.path.join(cfg.PATHS['climate_file'],
                                                    '**/*' + var + '_WFDE5_CRU+GPCC_' + '*_v2.0.nc')))
    else:
        paths_files = sorted(glob.glob(os.path.join(cfg.PATHS['climate_file'],
                                                '**/*' + var + '_WFDE5_CRU_' + '*_v2.0.nc')))

    files = []
    for path in paths_files:
        match = re.findall(r'\d+', path)
        if match:
            number = int(match[1])
            if int(y0) <= number <= int(y1):
                files.append(path)

    assert os.path.isfile(files[0])
    assert os.path.isfile(files[-1])
    print("File for beguinning and end of the time series exist")
    print(files[0])
    print(files[-1])
    return files


def _preprocess(x, i, j):
    """
    Preprocess to pass to xarray.open_mfdataset()
    :param x: Xarray Dataset
    :param i: index  for longitude
    :param j: index for latitude
    :return: a selected xarray Dataset for two specific coordinates
    """
    return x.isel(lat=j, lon=i)


def xropen_mfdataset(files,
                     lon=None,
                     lat=None):
    """
    Wrapper around xr.open_mfdataset() to pass a specific set of paths
    per climate variable
    :param files: files for each climate variable
    :param lon: longitude
    :param lat: latitude
    :return: a data set per variable concatenated by time and crop to
    the centre lat and lon of a glacier.
    """
    partial_func = partial(_preprocess, i=lon, j=lat)

    ds = xr.open_mfdataset(files,
                           concat_dim='time',
                           preprocess=partial_func,
                           engine='netcdf4',
                           combine='nested')

    return ds.load()


@entity_task(log, writes=['climate_historical_fsm'])
def process_wfde5_data(gdir,
                       y0=None,
                       y1=None):
    """
    Process WFDE5 meteorological variables and put them in an FSM ready format
    per glacier directory
    """
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=300).start()

    output_file_name = gdir.get_filepath('climate_historical_fsm')
    if os.path.exists(output_file_name):
        os.remove(output_file_name)

    lat = gdir.cenlat
    lon = gdir.cenlon

    if y0 is None:
        y0 = '1980'
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

    lwdown = 'LWdown'
    lwdown_fpaths = find_files_per_var(lwdown, y0=y0, y1=y1)
    psurf = 'PSurf'
    psurf_fpaths = find_files_per_var(psurf, y0=y0, y1=y1)
    qair = 'Qair'
    qair_fpaths = find_files_per_var(qair, y0=y0, y1=y1)
    rainf = 'Rainf'
    rainf_fpaths = find_files_per_var(rainf, y0=y0, y1=y1)
    snowf = 'Snowf'
    snowf_fpaths = find_files_per_var(snowf, y0=y0, y1=y1)
    swdown = 'SWdown'
    swdown_fpaths = find_files_per_var(swdown, y0=y0, y1=y1)
    tair = 'Tair'
    tair_fpaths = find_files_per_var(tair, y0=y0, y1=y1)
    wind = 'Wind'
    wind_fpaths = find_files_per_var(wind, y0=y0, y1=y1)

    # Prepare the data for Pool of workers
    paths = [lwdown_fpaths, psurf_fpaths, qair_fpaths, rainf_fpaths, snowf_fpaths, swdown_fpaths, tair_fpaths,
             wind_fpaths]
    ii = np.concatenate([[i]] * 8, axis=0)
    jj = np.concatenate([[j]] * 8, axis=0)

    d0 = datetime(int(y0), 1, 1, 0, 0, 0)
    d1 = datetime(int(y1), 12, 31, 23, 0, 0)
    delta = d1 - d0
    dimensions = (delta.days + 1) * 24

    coords = dict(time=(range(dimensions)), lon=None, lat=None)
    dlw = xr.DataArray(None, coords=coords, dims=("time",), name=lwdown, attrs=None)
    dsurf = xr.DataArray(None, coords=coords, dims=("time",), name=psurf, attrs=None)
    dqair = xr.DataArray(None, coords=coords, dims=("time",), name=qair, attrs=None)
    drainf = xr.DataArray(None, coords=coords, dims=("time",), name=rainf, attrs=None)
    dsnowf = xr.DataArray(None, coords=coords, dims=("time",), name=snowf, attrs=None)
    dswdown = xr.DataArray(None, coords=coords, dims=("time",), name=swdown, attrs=None)
    dtair = xr.DataArray(None, coords=coords, dims=("time",), name=tair, attrs=None)
    dwind = xr.DataArray(None, coords=coords, dims=("time",), name=wind, attrs=None)

    # Only 8 nodes at the time per glacier
    workers = 8
    with multiprocessing.Pool(processes=workers) as pool:
        dlw, dsurf, dqair, drainf, dsnowf, dswdown, dtair, dwind = pool.starmap(xropen_mfdataset,
                                                                                zip(paths,
                                                                                    ii,
                                                                                    jj)
                                                                               )
        pool.close()
        pool.join()

    # Merge all variables into a single data frame
    ds = xr.merge([dlw, dsurf, dqair, drainf, dsnowf, dswdown, dtair, dwind])

    ds.attrs = {'author': 'Beatriz Recinos and Richard Essery',
                'author_info': 'Big thaw OGGM-FSM',
                'ref_hgt': ref_hgt,
                'ref_pix_lat': ref_pix_lat,
                'ref_pix_lon': ref_pix_lon,
                'climate_source': 'WFDE5',
                'yr_0': y0,
                'yr_1': y1}

    ds.load().to_netcdf(gdir.get_filepath('climate_historical_fsm'),
                        mode='w',
                        format='NETCDF4',
                        engine='netcdf4')

    pbar.finish()


class FactorialSnowpackModel(MassBalanceModel):
    def __init__(self,
                 gdir,
                 filename='climate_historical_fsm',
                 input_filesuffix='',
                 mb=0.,
                 zmin=None,
                 zmax=None,
                 Nbnd=None):
        super(FactorialSnowpackModel, self).__init__()
        self.hemisphere = 'nh'
        self.valid_bounds = [-2e4, 2e4]  # in m

        # FSM layers
        self.zmin = zmin  # Centre of lowest elevation band (m)
        self.zmax = zmax  # Centre of highest elevation band (m)
        self.Nbnd = Nbnd  # Number of elevation bands
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
            dates = netCDF4.num2date(self.time, units=nc['time'].units, 
                    calendar=nc['time'].calendar)
            self.years = np.array([date.year for date in dates])
            self.months = np.array([date.month for date in dates])
        self._mb = mb


    def get_annual_mb(self, heights=None, year=None, fls=None, fl_id=None):

        if fls is None:
            raise RuntimeError(f'FSM requires flow band detail')
        else:
            areas = fls[0].bin_area_m2
            if heights is None:
                heights = fls[0].surface_h

        if year is not None:
            inds = np.where(self.years==year)
        else:
            inds = np.where(self.years > -99999)

        LW=self.LW[inds]
        Ps=self.Ps[inds]
        Qa=self.Qa[inds]
        Rf=self.Rf[inds]
        Sf=self.Sf[inds]
        SW=self.SW[inds]
        Ta=self.Ta[inds]
        Ua=self.Ua[inds]
        time=self.time[inds]
        Ntim=int(len(time))
        Nseg=int(len(areas))

        mb = FSM.fsmpy(self.Dice, self.Dmin, self.dz, LW, Ps,
                       Qa, Rf, Sf, SW, Ta, Ua,
                       areas, heights, self.albs, self.Dsnw, self.Nsnw, self.Sice,
                       self.Sliq, self.Tice, self.Tsnw, self.Tsrf, self.Nbnd,
                       self.Nice, self.Nsmx, Ntim, Nseg)

        # output is in kg / m^2 -- need to convert to m/s over a suitable baseline
        
        if year is None:
            baseline_y = np.max(self.years) - np.min(self.years) + 1
        else:
            baseline_y = 1

        rho = self.rho = cfg.PARAMS['ice_density']
        mb = (mb / baseline_y) / SEC_IN_YEAR / rho

        if min(heights) < self.zmin or max(heights) > self.zmax:
            raise RuntimeError(f'The heights provided are outside of the '
                                   f'elevation bands for FSM')

        else:

            func = interp1d(self.zbnd, mb)
            return func (heights)


    def is_year_valid(self, year):
        return self.ys <= year <= self.ye
