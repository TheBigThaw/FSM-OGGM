import numpy as np
import geopandas as gpd
from oggm import cfg, utils
from oggm import workflow, tasks
from FSM_oggm_MB import FactorialSnowpackModel, process_wfde5_data

cfg.initialize()

cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['mp_processes'] = 2
cfg.PARAMS['border'] = 80
cfg.initialize(logging_level='DEBUG')

reset=False
print('Reset is set to ', reset)
print('**Important set this to False to avoid '
      'resetting the glacier directory everytime this is ran!**')

cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-FSM-test', reset=reset)
#cfg.PATHS['working_dir'] = '/home/dgoldber/ice_models/oggm'
print('we are working here', cfg.PATHS['working_dir'])
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

fr = utils.get_rgi_region_file(11, version='62', reset=False)
gdf = gpd.read_file(fr)

## Selecting the glaciers that belong to the Rofental catchment
minlat = 46.690
maxlat = 47.170
minlon = 10.6
maxlon = 11.3
zmin=2050
zmax=3739
Nbnd=15
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

elevation_band_task_list = [
    tasks.simple_glacier_masks,
    tasks.elevation_band_flowline,
    tasks.fixed_dx_elevation_band_flowline,
    tasks.compute_downstream_line,
    tasks.compute_downstream_bedshape,
    tasks.gridded_attributes,
    tasks.gridded_mb_attributes,
]

for task in elevation_band_task_list:
    workflow.execute_entity_task(task, gdirs)

cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/boreal/WFDE5/'
cfg.PARAMS['baseline_climate'] = 'CUSTOM'

# placeholder until Dan fixes the climate preprocessing data
# for now we just copy and paste a file
#workflow.execute_entity_task(process_wfde5_data, gdirs, y0='1980', y1='2019')
#print ("DONE PROCESSING wfde5 data")

gdir = gdirs[0]

mass_balance = FactorialSnowpackModel(gdir, filename='climate_historical_fsm', zmin=zmin, zmax=zmax, Nbnd=Nbnd)
fls = gdir.read_pickle('inversion_flowlines')
mass_balance.get_annual_mb(fls=fls)

# Placeholder for next inversion steps until oggm changes
# apparent_mb_from_any_mb() task
tasks.apparent_mb_from_any_mb(gdir, mb_model=mass_balance, mb_years=np.unique(mass_balance.years))

workflow.calibrate_inversion_from_consensus(
    gdirs,
    apply_fs_on_mismatch=True,
    error_on_mismatch=True,  # if you're running many glaciers some might not work
    filter_inversion_output=True,  # this partly filters the over deepening due to
    # the equilibrium assumption for retreating glaciers (see. Figure 5 of Maussion et al. 2019)
    volume_m3_reference=None,  # here you could provide your own total volume estimate in m3
)

# finally create the dynamic flowlines
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)
print('all worked')




# mb_ts, zbnd = mbmod.get_annual_mb()
#import matplotlib.pyplot as plt
#mb = mb_ts / 40
#plt.plot(zbnd, mb, 'k')
#plt.xlim(2000, 4000)
#plt.xlabel('Elevation (m)')
#plt.ylabel('Annual mass balance (mm w.e.)')
#plt.savefig(os.path.join(cfg.PATHS['working_dir'], 'test.png'))
