#!/bin/bash

set -e

function build_pybind {

  git clone git@github.com:pybind/pybind11.git env/src/pybind11
  pushd env/src/pybind11

  cmake . -DCMAKE_INSTALL_PREFIX=../..
  make
  make install
  python3 setup.py install

  popd

}

function build_multem {

  git clone git@github.com:jmp1985/MULTEM.git env/src/MULTEM

}

# Setup developement environment
python3 -m venv env
source env/bin/activate

# Install the requirements
pip install -r requirements.txt

# Make a source directory
mkdir -p env/src

# Build dependancies
build_pybind
build_multem
