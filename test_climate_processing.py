import geopandas as gpd
from oggm import cfg, utils
from oggm import workflow, tasks
from IPython import embed

cfg.initialize()

cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['mp_processes'] = 2
cfg.PARAMS['border'] = 80

reset=False

cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-climate', reset=reset)
cfg.PARAMS['continue_on_error'] = True
cfg.PARAMS['use_compression'] = True
cfg.PARAMS['use_tar_shapefiles'] = True
cfg.PATHS['rgi_version'] = '62'
cfg.PARAMS['use_temp_bias_from_file'] = True
cfg.PARAMS['compress_climate_netcdf'] = False
cfg.PARAMS['store_model_geometry'] = True
cfg.PARAMS['store_fl_diagnostics'] = True

base_url = ('https://cluster.klima.uni-bremen.de/~oggm/'
            'gdirs/oggm_v1.6/L3-L5_files/2023.1/elev_bands/W5E5_w_data/')

fr = utils.get_rgi_region_file(11, version='62', reset=reset)
gdf = gpd.read_file(fr)

## Selecting the glaciers that belong to the Rofental catchment
minlat = 46.690
maxlat = 47.170
minlon = 10.6
maxlon = 11.3
rof = gdf[gdf['CenLat'].between(minlat, maxlat) & gdf['CenLon'].between(minlon, maxlon)]

rof = rof.sort_values('Area', ascending=False)

selection = rof[rof.Name == 'Hintereisferner']

if reset:
    gdirs = workflow.init_glacier_directories(selection,
                                          from_prepro_level=3,
                                          prepro_base_url=base_url,
                                          reset=reset)
else:
    gdirs = workflow.init_glacier_directories(selection)

cfg.PATHS['climate_file'] = '/home/brecinos/FSM/ASurf_WFDE5_CRU_v2.0.nc'
cfg.PARAMS['baseline_climate'] = 'CUSTOM'

embed()