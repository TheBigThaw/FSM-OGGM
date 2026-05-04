rm -r FSM.egg-info/  build/
rm fsm.pyf
rm *.so *.mod
python -m numpy.f2py FSM.f90 -m FSM -h fsm.pyf
FC="gfortran" python -m numpy.f2py -c fsm.pyf FSM.f90 --backend meson
#python setup.py build
#python setup.py install
