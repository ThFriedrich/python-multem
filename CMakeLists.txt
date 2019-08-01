cmake_minimum_required(VERSION 3.11.0)

# Set the project name
project(multem CXX CUDA)

# Need pybind11 for Python C/C++ extensions
find_package(pybind11 REQUIRED)
find_package(CUDA REQUIRED)

# Turn off LTO (has problems with cuda)
set(HAS_FLTO False)
set(HAS_FLTO_THIN False)

# Add a C/C++ extension
pybind11_add_module(multem_ext 
  src/multem/multem.cc
  src/multem/multem.cu)

# Ensure we are using C++11
target_compile_features(multem_ext PUBLIC cxx_std_11)

# Set the coverage compile and link flags for gcc or clang
set(COVERAGE_COMPILE_FLAGS -coverage)
set(COVERAGE_LINK_FLAGS    -coverage)

# Set the coverage compile and link flags in the debug build
target_compile_options(multem_ext 
  PUBLIC "$<$<CONFIG:DEBUG>:${COVERAGE_COMPILE_FLAGS}>")
target_link_options(multem_ext 
  PUBLIC "$<$<CONFIG:DEBUG>:${COVERAGE_LINK_FLAGS}>")

# Install the python extension
install(TARGETS multem_ext LIBRARY DESTINATION src/multem)
