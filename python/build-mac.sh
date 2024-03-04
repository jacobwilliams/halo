#
# Create a conda environment and build qr_mumps (serial version) on a mac
#

set -e    # quit on any errors
shopt -s expand_aliases

export CONDAENVDIR=$PWD/env
export QRMUMPSREPO=https://gitlab.com/qr_mumps/qr_mumps.git
export QRMUMPSSHA=9ea78d6630aea7a29e00c230d076b8bee0f48a1f
export QRMUMPSDIR=$PWD/qr_mumps
export INSTALLDIR=$QRMUMPSDIR/install

rm -rf $QRMUMPSDIR
rm -rf $CONDAENVDIR

# conda garbage:
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh    # this seems to work

conda env create -f ./environment.yml --prefix $CONDAENVDIR

conda activate $CONDAENVDIR

# pull the qr_mumps repo
git clone $QRMUMPSREPO
cd $QRMUMPSDIR
git reset --hard $QRMUMPSSHA

mkdir build
cd build

if [[ "$(uname)" == "Darwin" ]]; then
    cmake -Wno-dev -DARITH=d -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=ON -DQRM_ORDERING_METIS=ON  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DLAPACK_LIB_DIR=$CONDAENVDIR/lib -DLAPACK_LIB="libblas.dylib;liblapack.dylib" ..
else
    cmake -Wno-dev -DARITH=d -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=ON -DQRM_ORDERING_METIS=ON  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DLAPACK_LIB_DIR=$CONDAENVDIR/lib ..
fi

cmake --build . --parallel --config Release
cmake --install .
