#!/bin/bash
#Usage: make target CORE=[core] [options]
#Example targets:
#    ifort
#    gfortran
#    xlf
#    pgi
#    intel-xd2000 :: cptec/inpe ian xd2000
#Availabe Cores:
#    atmosphere
#    init_atmosphere
#    landice
#    ocean
#    seaice
#    sw
#    test
#Available Options:
#    DEBUG=true    - builds debug version. Default is optimized version.
#    USE_PAPI=true - builds version using PAPI for timers. Default is off.
#    TAU=true      - builds version using TAU hooks for profiling. Default is off.
#    AUTOCLEAN=true    - forces a clean of infrastructure prior to build new core.
#    GEN_F90=true  - Generates intermediate .f90 files through CPP, and builds with them.
#    TIMER_LIB=opt - Selects the timer library interface to be used for profiling the model. Options are:
#                    TIMER_LIB=native - Uses native built-in timers in MPAS
#                    TIMER_LIB=gptl - Uses gptl for the timer interface instead of the native interface
#                    TIMER_LIB=tau - Uses TAU for the timer interface instead of the native interface
#    OPENMP=true   - builds and links with OpenMP flags. Default is to not use OpenMP.
#    OPENACC=true  - builds and links with OpenACC flags. Default is to not use OpenACC.
#    USE_PIO2=true - links with the PIO 2 library. Default is to use the PIO 1.x library.
#    PRECISION=single - builds with default single-precision real kind. Default is to use double-precision.
#    SHAREDLIB=true - generate position-independent code suitable for use in a shared library. Default is false.

#module purge
#module load xpmem/0.2.119-1.3_gef379be13330
#module load PrgEnv-gnu/8.6.0
#module load craype-x86-turin
#module load cray-hdf5/1.14.3.3
#module load cray-netcdf/4.9.0.15
#module load cray-parallel-netcdf/1.12.3.15
#export NETCDF=/opt/cray/pe/netcdf/4.9.0.15/gnu/12.3
#export PNETCDF=/opt/cray/pe/parallel-netcdf/1.12.3.15/gnu/12.3

# Load modules:
module purge
module load PrgEnv-intel
module load craype-x86-turin
module load cray-hdf5/1.14.3.3
module load cray-netcdf/4.9.0.15
module load cray-parallel-netcdf/1.12.3.15
module load grads/2.2.1
module load cdo/2.4.2
module load METIS/5.1.0
module load cray-pals
module list

# Others variables:
export OMP_NUM_THREADS=1
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"
export PMIX_MCA_gds=hash
export MPI_PARAMS="-iface ib0 -bind-to core -map-by core"

export NETCDF=/opt/cray/pe/netcdf/4.9.0.15/INTEL/2023.2
export PNETCDF=/opt/cray/pe/parallel-netcdf/1.12.3.15/INTEL/2023.2
# PIO is not necessary for version 8.* If PIO is empty, MPAS Will use SMIOL
export PIO=
export LD_LIBRARY_PATH=$NETCDF/lib:$PNETCDF/lib:$PIO/lib:$LD_LIBRARY_PATH

DATE_TIME_NOW=$(date +"%Y%m%d%H%M")
MAKE_OUT_FILE="make_${DATE_TIME_NOW}.out"

#make clean CORE=init_atmosphere
#make -j 8 gfortran CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee -a ${MAKE_OUT_FILE}
#make -j 8 intel CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}
#make -j 8 intel2-xd2000 CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}
#make -j 8 intel-xd2000 CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single OPTIMIZATION_LEVEL=O1 FFLAGS_OPT=-O1 CFLAGS_OPT=-O1 CXXFLAGS_OPT=-O1 2>&1 | tee ${MAKE_OUT_FILE}

#make clean CORE=atmosphere
#make -j 8 gfortran CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}
#make -j 8 intel CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}
make -j 1 intel-xd2000 CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}

if [ -s "./init_atmosphere_model" ] && [ -e "./atmosphere_model" ]; then
    echo ""
    echo -e "\033[1;32m==>\033[0m Files init_atmosphere_model and atmosphere_model generated Successfully!"
    echo
else
    echo -e "\033[1;31m==>\033[0m !!! An error occurred during build. Check output"
    exit -1
fi

#-----
/bin/ln -fs ../tables/MP_THOMPSON_QIautQS_DATA.DBL   .
/bin/ln -fs ../tables/MP_THOMPSON_QRacrQG_DATA.DBL   .
/bin/ln -fs ../tables/MP_THOMPSON_QRacrQS_DATA.DBL   .
/bin/ln -fs ../tables/MP_THOMPSON_freezeH2O_DATA.DBL .
