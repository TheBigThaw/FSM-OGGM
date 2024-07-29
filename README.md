# FSM-OGGM coupling

This is a "draft" repository to build a mass balance class, which couples [The Factorial Snowpack Model](https://gmd.copernicus.org/articles/8/3867/2015/gmd-8-3867-2015.html) with the [Open Global Glacier Model](https://docs.oggm.org/en/latest/) and test the setup on Alpine Glaciers.


## How to install?

1. First clone the repository

``
git clone https://github.com/bearecinos/FSM-OGGM.git
``

2. Build a python environment compatible with `f2py`. For this, install an environment according to the following yml.

```
name: oggm_fsm
channels:
  - conda-forge
dependencies:
  - python==3.9.2
  - numpy==1.26.4
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
  - pip 
  - pip:
    - tables==3.9.2
    - joblib
    - progressbar2
```

A copy of this `environment_oggm_fsm.yml` file can be found under the FSM folder. Built your environment with mamba via:

```
mamba env create -f environment_oggm_fsm.yml
```

> At the moment we are fixing our numpy and python version due to an issue with [numpy.distutils.core and f2py](https://numpy.org/doc/stable/reference/distutils_status_migration.html#status-of-numpy-distutils-and-migration-advice). However, I think we can fix this by finding another way to write `setup.py` file and build the Fortran executable via a different method; I ran out of time to test another way of doing this. For more context see [f2py docs](https://numpy.org/doc/stable/f2py/index.html#f2py-user-guide-and-reference-manual) and this [blog](https://hackmd.io/@python-fortran-interface/SJ8kiUctd#f2py). Probably there is a better way!
 pytables is fixed too, so there are no conflicts with that numpy version and oggm test pass.

For the OGGM installation, I recommend to install it via Github and from my branch, as some test are broken in the main OGGM repo:
```
git clone https://github.com/bearecinos/oggm.git
cd oggm
pip install -e .
```

3. Built FSM python module with `f2py`.

```
cd FSM
```

**Dont forget to activate your oggm_fsm env and be inside the FSM folder!**. Then run the following commands:

```
python -m numpy.f2py FSM.f90 -m FSM -h fsm.pyf
python -m numpy.f2py -c fsm.pyf FSM.f90
python setup.py build
python setup.py install
```

4. Run OGGM test, ideally there should be one failing which can be ignored, see the following [issue](https://github.com/OGGM/oggm/issues/1714). Then test if FSM library has been installed correctly by opening python and importing the module.

```
import FSM
```

> Note: this repository code and documentation is a work in progress and might change alot dure to offline FSM development