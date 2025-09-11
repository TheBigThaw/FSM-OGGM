# FSM-OGGM coupling

This is a repository to build a mass balance class, which couples [The Factorial Snowpack Model](https://gmd.copernicus.org/articles/8/3867/2015/gmd-8-3867-2015.html) with the [Open Global Glacier Model](https://docs.oggm.org/en/latest/) and test the setup on Alpine Glaciers.

A subdirectory contains an experiment focused on the Roftental catchment, Austria, and scripts are provided for postprocessing.

## How to install?

1. First clone the repository

```
git clone https://github.com/bearecinos/FSM-OGGM.git
```

2. Build a python environment compatible with `f2py`. For this, install an environment according to the following yml. `meson` and `ninja` are important to compile the fortran executable.

```
name: oggm_fsm
channels:
  - conda-forge
dependencies:
dependencies:
  - python==3.11.0
  - numpy
  - jupyter
  - jupyterlab
  - scipy
  - pandas
  - shapely
  - matplotlib
  - Pillow
  - netcdf4
  - scikit-image
  - scikit-learn
  - configobj
  - xarray
  - pytest
  - dask
  - bottleneck
  - pyproj
  - cartopy
  - geopandas
  - rasterio
  - rioxarray
  - seaborn
  - salem
  - motionless
  - ipython
  - f90nml
  - pip
  - pip:
    - meson
    - ninja
    - tables
    - joblib
    - progressbar2
```

A copy of this lives in the `environment_oggm_fsm.yml` file, which can be found under the `FSM-OGGM/FSM` folder. Built your environment with mamba via:

```
mamba env create -f environment_oggm_fsm.yml
```

3. For the OGGM installation, I recommend to install it via Github:
```
git clone https://github.com/OGGM/oggm.git
cd oggm
pip install -e .
```

4. Don't foget to test your OGGM environment via:

```
mamba activate oggm_fsm
pytest.oggm  --disable-warnings
```

5. If all is well you are ready to build the FSM python module with, `meson`, `ninja` and `f2py`. You will also require a gfortran compiler (other versions such as intel fortran may work but have not been tested).

**Dont forget to activate your oggm_fsm env and be inside FSM-OGGM/FSM folder**. Then run `bash compil.sh`.

6. Then test if FSM-OGGM library has been installed correctly by opening python and importing the FSM module.
```
cd FSM-OGGM
python
```
Then in python do:
```
import FSM
```


## Rofental Experiment

The folder `slurm_scripts_rof` contains a number of python and shell scripts designed to implement an experiment for the Rofental catchment. Their use is explained here.

The `FactorialSnowpackModel` class is not "installable" -- the python source path must point to it. Before running these scripts, navigate to this directory and call the commands:

```
SCRIPT_DIR="$(pwd)"
ROOT_DIR="$(realpath "$SCRIPT_DIR/..")"
export PYTHONPATH="$ROOT_DIR:$PYTHONPATH"
```

There is a file `params.ini`, containing parameters for this experiment; we show it here:

```python
[General]
working_dir = None
reset = True

[OGGM]
use_multiprocessing = True
mp_processes = 20

[InputData]
climate_file = /exports/geos.ed.ac.uk/iceocean/WFDE5_rof/
glacier_rgi_id = None
y0 = 1980
y1 = 2019
catchment_path = /exports/geos.ed.ac.uk/iceocean/TBT/rofental_catchment_shp/02ee24e5-d25e-4aa6-b6a2-9d08a12992df_boundaries/deims_sites_boundariesPolygon.shp

[Output]
simulation_name = _climate_historical_fsm

[FSM_OGGM]
FSM_save_runoff = True
FSM_runoff_frequency = D
FSM_spinup = True
FSM_interpolate_bnds = False
FSM_Nbnds = 15
FSM_param_asmx = 0.85
FSM_param_asmn = 0.6
FSM_param_aice = 0.5
FSM_param_Plapse = 0.0
FSM_param_Pf = 1.2
FSM_param_Tlapse = 6.5e-3
FSM_param_sigmoidDscale = 0
```

It is necessary to set `working_dir` to an existing folder; this is where all output will be written. 

### Running the simulation with FSM

To run the simulation, call

`python test_fsm_rofental.py params.ini`

This should take several minutes to nearly half an hour depending on how many processors are used. If finished successfully, there will be outputs in `[workdir]/per_glacier`. Look for `fl_diagnostics_climate_historical_fsm.nc` and `gridded_simulation_climate_historical_fsm.nc` in one of the glacier directories, their presence indicates things went correctly.

### Glacier outlines

Annual shapefiles are then generated with

`python output_area_change_shapefiles.py params.ini`.

Shapefiles are saved in `[workdir]/area_evolution/[simulation_name].

### Glacier distributed thickness

Gridded maps of thickness are generated with 

`python output_distributed_thickness_and_runoff.py params.ini`. This script saves outputs (as netcdf files) in `[workdir]/all_simulations_merged[simulation_name].nc`. Note that this function also produces a *runoff file* in the folder `run_off_terminus_position`. This runoff file gives snow and ice melt runoff by *day* and by *RGI code*.

**IMPORTANT** this script will fail if you do not have a NASA EarthData credentials file. *To creat one*:

- First, create an account at [https://urs.earthdata.nasa.gov/] if you do not have one already.
- Next, call the command-line interface function `oggm_netrc_credentials` and follow prompts.
- If calling `oggm_netrc_credentials` does *not* work, you can call the python function directly.
-   Navigate to `oggm/oggm/cli/`
-   Open a python or Ipython terminal
-   Call `from netrc_credentials import cli` and then call `netrc_credentials.cli()`. Follow prompts.

### Glacier Terminus Positions

Finally, call `output_terminus_position_to_runoff_file.py params.ini`. This script calculated terminus position of each glacier at the *start* of each year. These positions (in lat/lon) are then saved into the `run_off_terminus_position.nc` file. While the runoff file is daily, position is found only at the start of each year, so only the entries in `run_off_terminus_position.nc` that fall on 1 Jan have nonempty values.

Note that `output_area_change_shapefiles.py` can be called directly after `test_fsm_rofental.py` or after `output_distributed_thickness_and_runoff.py` and/or `output_terminus_position_to_runoff_file.py`. `output_terminus_position_to_runoff_file.py` has dependency on `output_distributed_thickness_and_runoff.py`. All output scripts require `test_fsm_rofental.py` to be run first.
