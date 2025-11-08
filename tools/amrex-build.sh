#!/bin/bash
# This script builds AMReX

# Usage: amrex-build.sh -b <build_dir> -c <compiler> -i <install_dir>

# Set the build argument paased to the script to variable BUILDIR
# collecting the arguments passed to the script

while getopts "b:c:i:" opt; do
    case $opt in
        b) BUILDIR=$OPTARG;;
        c) COMPILER=$OPTARG;;
        i) INSTALL_DIR=$OPTARG;;
        *) echo "Invalid option"; exit 1;;
    esac
done

# set the correct CC and CXX compiler based on the compiler argument
export GPU_BACKEND=NONE
export AMReX_DIFFERENT_COMPILER=OFF
if [ "$COMPILER" == "gcc" ]; then
    export CC=gcc
    export CXX=g++
    export FC=gfortran
    export F77=gfortran

elif [ "$COMPILER" == "gcc-9" ]; then
    export CC=gcc-9
    export CXX=g++-9
    export FC=gfortran-9
    export F77=gfortran-9

elif [ "$COMPILER" == "clang" ]; then
    export CC=clang
    export CXX=clang++
    export FC=flang
    export F77=flang
elif [ "$COMPILER" == "intel" ]; then
    export CC=icx
    export CXX=icpx
    export FC=ifx
    export F77=ifx
    export MPI_CXX=mpiicpx
    export MPI_C=mpiicx
elif [ "$COMPILER" == "nvhpc" ]; then
    export CC=nvc
    export CXX=nvc++
    export FC=nvfortran
    export F77=nvfortran
    export GPU_BACKEND=CUDA
    export AMReX_DIFFERENT_COMPILER=ON
else
    echo "Invalid compiler specified. Please use gcc, clang, nvhpc or intel."
    exit 1
fi
# check if the build directory is set
if [ -z "$BUILDIR" ]; then
    echo "Build directory not set. Using default: build"
    BUILDIR="build"
fi

# check if the install directory is set
if [ -z "$INSTALL_DIR" ]; then
    echo "Install directory not set. Using default: install"
    INSTALL_DIR="install"
fi
# check if the compiler is set
if [ -z "$COMPILER" ]; then
    echo "Compiler not set. Using default: g++"
    COMPILER="g++"
fi

# Set the compiler argument to a variable COMPILER
#
AMReX_CONFIG=1
AMReX_Root=${PWD}
AMReX_INSTALL_PREFIX=${AMReX_Root}/${INSTALL_DIR}/${COMPILER}
AMReX_MAKE=1
AMReX_CLOBBER=0
AMReX_TEST=0
AMReX_INSTALL=1
AMReX_VERBOSE=0

ENABLE_MPI=1
ENABLE_OpenMP=0
#MPI_PREFIX=${OMPI}

# 1, 2, 3
SPACEDIM=3

# FLOAT, DOUBLE
PRECISION=DOUBLE

#debug, release
CMAKE_BUILD_TYPE=release

# cmake, ctest executables
CMAKE=cmake
CTEST=ctest

# Relative to AMReX_Root
Build_Dir=${BUILDIR}

# cores for build
AMReX_NPROCS=12

while getopts "ab:cfin:mtv" flag
do
  case $flag in
    a) AMReX_CONFIG=1; AMReX_MAKE=1; AMReX_TEST=1; AMReX_INSTALL=1;;
    b) Build_Dir=${OPTARG};;
    c) AMReX_CONFIG=1;;
    f) AMReX_CLOBBER=1;;
    i) AMReX_INSTALL=1;;
    n) AMReX_NPROCS=${OPTARG};;
    m) AMReX_MAKE=1;;
    t) AMReX_TEST=1;;
    v) AMReX_VERBOSE=1;;
  esac
done

AMReX_Build=${AMReX_Root}/${Build_Dir}

echo AMReX_Root=$AMReX_Root
echo AMReX_Build=$AMReX_Build
echo AMReX_CONFIG=$AMReX_CONFIG
echo AMReX_MAKE=$AMReX_MAKE
echo AMReX_TEST=$AMReX_TEST
echo AMReX_INSTALL=$AMReX_INSTALL
echo AMReX_NPROCS=$AMReX_NPROCS

if [ $AMReX_CLOBBER -eq 1 ]; then
    rm -rf ${AMReX_Build}
fi

if [ $AMReX_CONFIG -eq 1 ]; then
    mkdir -p ${AMReX_Build}
    cd ${AMReX_Build}

    CMAKE_ARGS="-D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} \
                -D ENABLE_Config_Report:BOOL=ON \
                -D BL_SPACEDIM:INT=${SPACEDIM} \
                -D BL_PRECISION:STRING=${PRECISION} \
                -D AMReX_ENABLE_TESTS:BOOL=ON \
                -D AMReX_AMRDATA:BOOL=ON\
                -D AMReX_AMRLEVEL:BOOL=ON\
                -D AMReX_BOUND_CHECK:BOOL=ON\
                -D AMReX_FORTRAN:BOOL=ON\
                -D AMReX_FORTRAN_INTERFACES:BOOL=ON\
                -D AMReX_BUILD_SHARED_LIBS:BOOL=ON\
                -D AMReX_EB:BOOL=ON\
                -D AMReX_GPU_BACKEND:STRING=${GPU_BACKEND} \
                -D AMReX_DIFFERENT_COMPILER:BOOL=${AMReX_DIFFERENT_COMPILER}\
                -D CMAKE_INSTALL_PREFIX:FILEPATH=${AMReX_INSTALL_PREFIX} "

    if [ ${AMReX_VERBOSE} -eq 1 ]; then
        CMAKE_ARGS=${CMAKE_ARGS} " --debug-output"
    fi

    if [ ${ENABLE_MPI} -eq 1 ]; then
        CMAKE_ARGS="${CMAKE_ARGS} -D ENABLE_MPI:BOOL=${ENABLE_MPI} \
                                   -D MPI_PREFIX:FILEPATH=${MPI_PREFIX}"

        if [ ${AMReX_TEST} -eq 1 ]; then
            CMAKE_ARGS="${CMAKE_ARGS} -D MPI_EXEC:FILEPATH=${MPI_PREFIX}/bin/mpiexec \
                                      -D MPI_EXEC_NUMPROCS_FLAG:STRING=-np \
                                      -D MPI_EXEC_ARGS:STRING='-mca mpi_yield_when_idle 1'"
        fi
    fi

    if [ ${ENABLE_OpenMP} -eq 1 ]; then
        CMAKE_ARGS="${CMAKE_ARGS} -D ENABLE_OpenMP:BOOL=${ENABLE_OpenMP}"
    fi

    ${CMAKE} ${CMAKE_ARGS} ${AMReX_Root}

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

if [ $AMReX_MAKE -eq 1 ]; then

    cd ${AMReX_Build}

    MAKE_ARGS=

    if [ ${AMReX_NPROCS} -ne 1 ]; then
        MAKE_ARGS="${MAKE_ARGS} -j ${AMReX_NPROCS}"
    fi

    if [ ${AMReX_VERBOSE} -eq 1 ]; then
        MAKE_ARGS="${MAKE_ARGS} VERBOSE=ON"
    fi

    make ${MAKE_ARGS}

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

if [ $AMReX_INSTALL -eq 1 ]; then
    cd ${AMReX_Build}

    make install

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

if [ $AMReX_TEST -eq 1 ]; then
    cd ${AMReX_Build}

    ${CTEST} --timeout 60 --output-on-failure

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

