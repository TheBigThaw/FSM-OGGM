# FSM-OGGM coupling

This is a "draft" repository to build a mass balance class, which couples [The Factorial Snowpack Model](https://gmd.copernicus.org/articles/8/3867/2015/gmd-8-3867-2015.html) with the [Open Global Glacier Model](https://docs.oggm.org/en/latest/) and test the setup on Alpine Glaciers.


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


7. Now you can run `test_fsm_oggm.py` via:
```
cd FSM-OGGM
python test_fsm_oggm.py
```

**Important**: once you ran the test once, make sure to set reset=False, so you don't have to produce again and again the glacier directory and download the data for the RGI, pre-process glacier dirs etc...

> Note: this repository code and documentation is a work in progress and might change alot dure to offline FSM development
