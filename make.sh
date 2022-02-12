#!/bin/bash

IDIR=${HOME}/build/include/coin-or
LDIR=${HOME}/build/lib
MEXBUILDINGBLOCKS=/mnt/c/repos/MexBuildingBlocks/

MATLABROOT=/opt/MATLAB/R2020b
MEX=${MATLABROOT}/bin/mex

${MEX} -output ipopt -v \
        -I${IDIR} \
        -L${LDIR} \
        -I${MEXBUILDINGBLOCKS} \
        CXXFLAGS='-std=c++17 -fPIC -fexceptions -fno-omit-frame-pointer -pthread' \
        LDFLAGS="-Wl,-rpath -Wl,${LDIR} -Wl,-rpath -Wl,${HOME}/OpenBLAS-static-libgcc -static-libstdc++" \
        -l:libipopt.a \
        -l:libcoinhsl.a \
        -L${HOME}/OpenBLAS -l:libopenblas.a -lgfortran \
        -L/usr/lib/x86_64-linux-gnu -lmetis \
        main.cpp

${MATLABROOT}/bin/matlab -nodisplay -nosplash -nojvm -r "hs71"