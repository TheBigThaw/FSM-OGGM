import argparse
import sys
import geopandas as gpd
import xarray as xr
import configparser
from oggm import cfg, utils
from oggm import workflow, tasks
from oggm.sandbox import distribute_2d
from FSM_oggm_MB import FactorialSnowpackModel, process_wfde5_data, fsm_flowline_model_run

def main(cfg_path):

    # ----------------------
    # 1) Read configuration file
    # ----------------------
    cp = configparser.ConfigParser()
    cp.read(cfg_path)

    gen  = cp['General']
    oggm = cp['OGGM']
    fsm  = cp['FSM_OGGM']
    inp  = cp['InputData']
    outp = cp['Output']

    # ----------------------
    # 2) Parse general settings
    # ----------------------
    working_dir = gen.get('working_dir')          # string, may be blank
    reset       = gen.getboolean('reset')         # bool

    # ----------------------
    # 3) Initialize OGGM core
    # ----------------------
    cfg.initialize(logging_level='DEBUG')
    cfg.PATHS['working_dir'] = utils.mkdir(working_dir, reset=reset)
    print(f"Working directory: {cfg.PATHS['working_dir']}")
    print(f"Reset = {reset}  (set False to preserve existing directories)")

    # ----------------------
    # 4) OGGM parameters
    # ----------------------
    cfg.PARAMS['use_multiprocessing'] = oggm.getboolean('use_multiprocessing')
    cfg.PARAMS['mp_processes']        = oggm.getint('mp_processes')
    cfg.PARAMS['border']              = oggm.getint('border', fallback=80)

    # ----------------------
    # 5) FSM parameters
    # ----------------------
    cfg.PARAMS['FSM_save_runoff']        = fsm.getboolean('FSM_save_runoff')
    cfg.PARAMS['FSM_runoff_frequency']   = fsm.get('FSM_runoff_frequency')
    cfg.PARAMS['FSM_spinup']             = fsm.getboolean('FSM_spinup')
    cfg.PARAMS['FSM_interpolate_bnds']   = fsm.getboolean('FSM_interpolate_bnds')
    # note: nbnds can only be set if interpolate_bnds is True
    cfg.PARAMS['FSM_Nbnds']              = fsm.getint('FSM_Nbnds')
    cfg.PARAMS['FSM_param_asmx']         = fsm.getfloat('FSM_param_asmx')
    cfg.PARAMS['FSM_param_asmn']         = fsm.getfloat('FSM_param_asmn')
    cfg.PARAMS['FSM_param_aice']         = fsm.getfloat('FSM_param_aice')
    cfg.PARAMS['FSM_param_Plapse']       = fsm.getfloat('FSM_param_Plapse')
    cfg.PARAMS['FSM_param_Pf']           = fsm.getfloat('FSM_param_Pf')
    cfg.PARAMS['FSM_param_Tlapse']       = fsm.getfloat('FSM_param_Tlapse')
    cfg.PARAMS['FSM_param_sigmoidDscale']= fsm.getint('FSM_param_sigmoidDscale')

    FactorialSnowpackModel.create_nml(reset=reset)

    # ----------------------
    # 6) Climate & I/O paths
    # ----------------------
    cfg.PATHS['climate_file']   = inp.get('climate_file')
    cfg.PARAMS['baseline_climate'] = 'CUSTOM'
    catchment_path = inp.get('catchment_path')

    # ----------------------
    # 7) OGGM run setup
    # ----------------------
    y0 = inp.getint('y0')
    y1 = inp.getint('y1')
    rgi_id = inp.get('glacier_rgi_id')
    simulation_name = inp.get('simulation_name')

    _doc = ('A netcdf file containing dates and ' +
            'ice-based and snow-based runoff volume ' +
            'for each date interval')
    cfg.BASENAMES['FSM_runoff'] = ('FSM_runoff.nc', _doc)

    print('Reset is set to ', reset)
    print('**Important set this to False to avoid '
          'resetting the glacier directory everytime this is ran!**')

    # this sets a temporary working directory. if you want to use a permanent
    # directory then uncomment and adapt the following line.
    cfg.PATHS['working_dir'] = utils.mkdir(working_dir, reset=reset)
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

    rof_shp = gpd.read_file(catchment_path)

    rof_sel = gdf.clip(rof_shp)
    rof_sel = rof_sel.sort_values('Area', ascending=False)

    if rgi_id in ('None', '', None):
        rgi_id = None

    if rgi_id:
        selection = rof_sel[rof_sel.RGIId == rgi_id]
    else:
        selection = rof_sel

    if reset:
        gdirs = workflow.init_glacier_directories(selection,
                                                  from_prepro_level=3,
                                                  prepro_base_url=base_url,
                                                  reset=reset,
                                                  force=reset)
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

    # Distribute
    workflow.execute_entity_task(tasks.distribute_thickness_per_altitude, gdirs)

    # Test that we have at least 21 variables on gridded_data.nc
    gdir = gdirs[0]
    with xr.open_dataset(gdir.get_filepath('gridded_data')) as ds:
        ds = ds.load()
        assert len(ds.count().variables.keys()) == 21


    workflow.execute_entity_task(process_wfde5_data, gdirs, y0=str(y0), y1=str(y1))
    print("DONE PROCESSING wfde5 data")

    workflow.execute_entity_task(tasks.apparent_mb_from_any_mb,
                                 gdirs,
                                 mb_model_class=FactorialSnowpackModel)

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

    workflow.execute_entity_task(fsm_flowline_model_run, gdirs,
                                 climate_filename='climate_historical_fsm',
                                 output_filesuffix='_climate_historical_fsm',
                                 ys=y0, ye=y1)

    print("DONE running FSM")


    workflow.execute_entity_task(distribute_2d.add_smoothed_glacier_topo, gdirs)
    workflow.execute_entity_task(distribute_2d.assign_points_to_band, gdirs)
    workflow.execute_entity_task(distribute_2d.distribute_thickness_from_simulation,
                                 gdirs, input_filesuffix=simulation_name)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: test_fsm_rofental.py <config.ini>")
        sys.exit(1)
    main(sys.argv[1])
