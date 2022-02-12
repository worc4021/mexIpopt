#!/bin/bash

# FROM ubuntu:18.04

# USER root

# RUN apt update
# RUN apt install -y build-essential 
# RUN apt install -y gfortran 
# RUN apt install -y cmake
# RUN apt install -y patch 
# RUN apt install -y wget 
# RUN apt install -y pkg-config 
# RUN apt install -y liblapack-dev 
# RUN apt install -y libmetis-dev
# RUN apt install -y git
# RUN apt install -y gcc-8 g++-8
# RUN apt install -y file
# RUN apt install -y xorg
# RUN apt install -y unzip
# RUN apt install -y libselinux1
# RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 700 --slave /usr/bin/g++ g++ /usr/bin/g++-7
# RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 800 --slave /usr/bin/g++ g++ /usr/bin/g++-8

# WORKDIR /home

# WORKDIR /home
# this folder is mounted into /source

PWD=$(pwd)
cd ${HOME}

echo "Deleting build output folder"
rm -rf ${HOME}/build
mkdir ${HOME}/build

MATLABROOT=/opt/MATLAB/R2020b

# Not sure why this needs to be exported but simply having a local variable does not work.
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${HOME}build/lib/pkgconfig/



echo "Deleting repos"
rm -rf ${HOME}/ThirdParty-HSL
rm -rf ${HOME}/Ipopt
rm -rf ${HOME}/OpenBLAS

git clone https://github.com/xianyi/OpenBLAS.git
cd ${HOME}/OpenBLAS
make
cd ${HOME}

git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
# download coinhsl (https://www.hsl.rl.ac.uk/ipopt/) which is free for academia and put it in ${HOME}/coinhsl or change the line below.
cp --recursive ${HOME}/coinhsl ${HOME}/ThirdParty-HSL/coinhsl
cd ${HOME}/ThirdParty-HSL
./configure --prefix=${HOME}/build --with-metis-lflags="-L/usr/lib/x86_64-linux-gnu -lmetis" --with-metis-cflags=-I/usr/include --enable-static --with-lapack-lflags="-L${HOME}/OpenBLAS -lopenblas"
make -j $(nproc) install

cd ${HOME}
git clone https://github.com/coin-or/Ipopt.git
cd ${HOME}/Ipopt
git checkout origin/stable/3.14

./configure --prefix=${HOME}/build --with-hsl-cflags="$(pkg-config --cflags coinhsl)" --with-hsl-lflags="$(pkg-config --libs coinhsl)" --disable-linear-solver-loader --disable-sipopt --disable-java --enable-static --with-lapack-lflags="-L${HOME}/OpenBLAS -lopenblas"
make -j $(nproc) install


cd ${PWD}