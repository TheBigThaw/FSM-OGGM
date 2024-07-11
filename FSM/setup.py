from numpy.distutils.core import Extension, setup

FSM = Extension(name = 'FSM',
                  sources = ['FSM.f90'])
if __name__ == "__main__":
    setup(name = 'FSM', ext_modules = [ FSM ])
