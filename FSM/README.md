# Built FSM python module with f2py


```
python -m numpy.f2py FSM.f90 -m FSM -h fsm.pyf
python -m numpy.f2py -c fsm.pyf FSM.f90
python setup.py build
python setup.py install
```

test by opening python and importing the module

```
import FSM
```