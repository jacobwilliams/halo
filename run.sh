#!/bin/bash

#
# This is meant to be used after running ../python/build-mac.sh
# (if you are using the qr_mumps solver). It will activate
# the conda environment and compile/link with qr_mumps and run
# the example file. uses the conda-forge compilers (gfortran, etc).
#

set -e    # quit on any errors
shopt -s expand_aliases

# export CONDAENVDIR=$PWD/env

# building on mac

export PYTHON_DIR=$PWD/python
export CONDA_ENV_DIR=$PYTHON_DIR/env

# conda garbage:
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh    # this seems to work

conda activate $CONDA_ENV_DIR

export DYLD_LIBRARY_PATH=$PYTHON_DIR/qr_mumps/install/lib:$DYLD_LIBRARY_PATH

fpm run halo_solver --profile release --flag "-DWITH_QRMUMPS -fopenmp $PYTHON_DIR/qr_mumps/install/lib/libdqrm.dylib $PYTHON_DIR/qr_mumps/install/lib/libqrm_common.dylib -I$PYTHON_DIR/qr_mumps/install/include -rpath $PYTHON_DIR/qr_mumps/install/lib" -- examples/example_sparse.json
