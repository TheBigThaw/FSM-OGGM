# meson.build
project('FSM', 'fortran', default_options: ['fortran_std=legacy'])

# Define the Fortran sources
fortran_sources = files('FSM/FSM.f90')

# Set up Python dependency
python = import('python').find_installation('python3')

## TODO: this needs a c wrapper.

# Create the Fortran-based shared library for the Python module
fsm_module = python.extension_module(
    'FSM',               # Name of the Python module
    fortran_sources,     # Fortran source files
    install: true        # Optional: whether to install the module
)
