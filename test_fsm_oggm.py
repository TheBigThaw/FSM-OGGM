import numpy as np
import geopandas as gpd
import pandas as pd
from oggm import cfg, utils
from oggm import workflow, tasks
from FSM_oggm_MB import FactorialSnowpackModel, process_wfde5_data, fsm_flowline_model_run, process_metum_data

cfg.initialize(logging_level='DEBUG')


# if multiprocessing is set to True, then 
# pooling will not be used for climate data
cfg.PARAMS['use_multiprocessing'] = False
cfg.PARAMS['mp_processes'] = 24
cfg.PARAMS['border'] = 80
cfg.PARAMS['FSM_interpolate_bnds'] = False
cfg.PARAMS['FSM_save_runoff'] = True
cfg.PARAMS['FSM_runoff_frequency'] = 'D'
#cfg.PARAMS['FSM_param_asmx'] = .99
_doc = ('A netcdf file containing dates and ' + 
        'ice-based and snow-based runoff volume ' + 
        'for each date interval')
cfg.BASENAMES['FSM_runoff'] = ('FSM_runoff.nc',_doc)

# the following if True means FSM will be run with a fixed # of columns
# independent on the number of glacier sectoins/elev bands. These
# will be spaced evenly over the elev range (which will include the 
# downstream region - so im not sure it is a good idea to use at all). 
# The number of columns/bands can be set in 
# the FactorialSnowpackModel constructor as a kwarg, or through the 
# cfg.PARAMS['FSM_Nbnds'] parameter (kwarg overwrites global param)
# or has a default of 15
cfg.PARAMS['FSM_interpolate_bnds'] = False
cfg.PARAMS['FSM_Nbnds'] = 15

# if True, this will run FSM for one year when FactorialSnowpackModel is 
# initiated, and the results will be the saved "initial state"
cfg.PARAMS['FSM_spinup'] = True

# Here is how an FSM parameter (asmx) is set. this will create a 
# namelist entry with the value equal to the default
cfg.PARAMS['FSM_param_asmx'] = .85

# this is necessary to create a nlst file in the present directory, which is the
# one FSM will read. would be good to enable to pass a path to the namelist to FSM
FactorialSnowpackModel.create_nml(reset=True)

reset=True
print('Reset is set to ', reset)
print('**Important set this to False to avoid '
      'resetting the glacier directory everytime this is ran!**')

# this sets a temporary working directory. if you want to use a permanent
# directory then uncomment and adapt the following line.
cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-FSM-test', reset=reset)
#cfg.PATHS['working_dir'] = '/exports/geos.ed.ac.uk/iceocean/dgoldber/FSM-OGGM'
print('we are working here', cfg.PATHS['working_dir'])

# bespoke path -- needs to be reset
cfg.PATHS['climate_file'] = '/exports/geos.ed.ac.uk/iceocean/WFDE5_rof/'
cfg.PATHS['metum_climate_file'] = '/exports/geos.ed.ac.uk/iceocean/dgoldber/FSM-OGGM/metum_temp/'
cfg.PARAMS['baseline_climate'] = 'CUSTOM'

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
ds_rof = pd.read_csv('rof_ids',header=None)

# By default only Hintereisferner is modeled.
# To model all values in the attached list uncomment below.
#selection = ds_rof[0].values.tolist()


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


print('multiprocessing' + str(cfg.PARAMS['use_multiprocessing']))

for task in elevation_band_task_list:
    workflow.execute_entity_task(task, gdirs)


workflow.execute_entity_task(process_metum_data, gdirs, y0='1999', y1='2001')
print ("DONE PROCESSING metum data")
workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs, mb_model_class=FactorialSnowpackModel)

workflow.calibrate_inversion_from_consensus(
    gdirs,
    apply_fs_on_mismatch=True,
    error_on_mismatch=True,  # if you're running many glaciers some might not work
    filter_inversion_output=True,  # this partly filters the over deepening due to
#    # the equilibrium assumption for retreating glaciers (see. Figure 5 of Maussion et al. 2019)
    volume_m3_reference=None,  # here you could provide your own total volume estimate in m3
)

# finally create the dynamic flowlines
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)

workflow.execute_entity_task(fsm_flowline_model_run,gdirs,
                             climate_filename='climate_historical_fsm',
                             ys=1981, ye=2019)

print('all worked')


