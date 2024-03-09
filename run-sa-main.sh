#!/bin/bash

#
# Script for the halo_optimizer
#

set -e    # quit on any errors
shopt -s expand_aliases

export PYTHON_DIR=$PWD/python
export CONDA_ENV_DIR=$PYTHON_DIR/env

# conda garbage:
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh    # this seems to work
conda activate $CONDA_ENV_DIR

if [[ "$(uname)" == "Darwin" ]]; then
    # have to do this for qr_mumps:
    export DYLD_LIBRARY_PATH=$PYTHON_DIR/qr_mumps/install/lib:$DYLD_LIBRARY_PATH

    # compile and run:

    # with qr_mumps:
    fpm run halo_optimizer --profile release --flag "-DWITH_QRMUMPS -fopenmp $PYTHON_DIR/qr_mumps/install/lib/libdqrm.dylib $PYTHON_DIR/qr_mumps/install/lib/libqrm_common.dylib -I$PYTHON_DIR/qr_mumps/install/include -rpath $PYTHON_DIR/qr_mumps/install/lib" -- $HALO_CONFIG_FILE

    # without qr_mumps + real128 kinds:
    # fpm run halo_solver --profile release --flag "-DREAL128 -fopenmp" -- $HALO_CONFIG_FILE

else

    # have to do this for qr_mumps:
    export LD_LIBRARY_PATH=$PYTHON_DIR/qr_mumps/install/lib:$LD_LIBRARY_PATH

    # compile the run:
    fpm run halo_optimizer --profile release --flag "-DWITH_QRMUMPS -fopenmp $CONDA_ENV_DIR/lib/libmetis.so $PYTHON_DIR/qr_mumps/install/lib/libdqrm.so $PYTHON_DIR/qr_mumps/install/lib/libqrm_common.so -I$PYTHON_DIR/qr_mumps/install/include" -- $HALO_CONFIG_FILE

fi

