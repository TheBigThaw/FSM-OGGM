import geopandas as gpd
from oggm import cfg, utils
from oggm import workflow, tasks
from FSM_oggm_MB import FactorialSnowpackModel, process_wfde5_data
from IPython import embed
from oggm.core.flowline import FileModel

cfg.initialize()

cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['mp_processes'] = 2
cfg.PARAMS['border'] = 80

reset=True
print('Reset is set to ', reset)
print('**Important set this to False to avoid '
      'resetting the glacier directory everytime this is ran!**')

cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-FSM-test', reset=reset)
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

fr = utils.get_rgi_region_file(11, version='62', reset=reset)
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

# Tested tasks
task_list = [
    tasks.compute_downstream_line,
    tasks.compute_downstream_bedshape,
    tasks.gridded_attributes,
    tasks.gridded_mb_attributes,
]
for task in task_list:
    workflow.execute_entity_task(task, gdirs)

# Distribute
workflow.execute_entity_task(tasks.distribute_thickness_per_altitude, gdirs)

cfg.PATHS['climate_file'] = '/exports/csce/datastore/geos/groups/boreal/WFDE5/'
cfg.PARAMS['baseline_climate'] = 'CUSTOM'


if (__name__ == '__main__') and reset:
    workflow.execute_entity_task(process_wfde5_data, gdirs, y0='1980', y1='2019')

print ("DONE PROCESSING wfde5 data")

gdir = gdirs[0]

mass_balance = FactorialSnowpackModel(gdir, filename='climate_historical_fsm', zmin=zmin, zmax=zmax, Nbnd=Nbnd)
gdir.get_inversion_flowline_hw()
fls = gdir.read_pickle('model_flowlines')
print(mass_balance.get_annual_mb(year=1990,fls=fls))

print('all worked')




# mb_ts, zbnd = mbmod.get_annual_mb()


#import matplotlib.pyplot as plt
#mb = mb_ts / 40
#plt.plot(zbnd, mb, 'k')
#plt.xlim(2000, 4000)
#plt.xlabel('Elevation (m)')
#plt.ylabel('Annual mass balance (mm w.e.)')
#plt.savefig(os.path.join(cfg.PATHS['working_dir'], 'test.png'))
