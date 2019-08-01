cmake_minimum_required(VERSION 3.11.0)

# Set the project name
project(python-multem)

# Need pybind11 for Python C/C++ extensions
find_package(pybind11 REQUIRED)

# Add a C/C++ extension
pybind11_add_module(multem_ext src/multem/ext.cc)

# Set the coverage compile and link flags for gcc or clang
set(COVERAGE_COMPILE_FLAGS -coverage)
set(COVERAGE_LINK_FLAGS    -coverage)

# Set the coverage compile and link flags in the debug build
target_compile_options(multem_ext 
  PUBLIC "$<$<CONFIG:DEBUG>:${COVERAGE_COMPILE_FLAGS}>")
target_link_options(multem_ext 
  PUBLIC "$<$<CONFIG:DEBUG>:${COVERAGE_LINK_FLAGS}>")

# Install the python extension
install(TARGETS multem_ext LIBRARY DESTINATION multem)
