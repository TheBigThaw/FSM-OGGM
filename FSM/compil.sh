rm -r FSM.egg-info/ *.so fsm.pyf
python -m numpy.f2py FSM.f90 -m FSM -h fsm.pyf
FC="gfortran" python -m numpy.f2py -c fsm.pyf FSM.f90 --backend meson
python setup.py build
python setup.py install
