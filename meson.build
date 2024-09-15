# meson.build
project('FSM_python_module', 'fortran')

# Define the Fortran sources
fortran_sources = files('FSM/FSM.f90')

# Set up Python dependency
python = import('python').find_installation('python3')

# Create the Fortran-based shared library for the Python module
FSM = python.extension_module(
    'FSM',  # Name of the module
    fortran_sources,   # Fortran source files
    install: true      # Optional: whether to install the module
)
