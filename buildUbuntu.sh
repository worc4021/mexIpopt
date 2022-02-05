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

cd /home

echo "Deleting build output folder"
rm -rf /home/build
mkdir /home/build

# Not sure why this needs to be exported but simply having a local variable does not work.
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/home/build/lib/pkgconfig/

echo "Deleting repos"
rm -rf /home/ThirdParty-HSL
rm -rf /home/Ipopt

git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
cp --recursive /source/coinhsl-archive-2014.01.17/ /home/ThirdParty-HSL/coinhsl
cd /home/ThirdParty-HSL
./configure --prefix=/home/build --with-metis-lflags="-L/usr/lib/x86_64-linux-gnu -lmetis" --with-metis-cflags=-I/usr/include --enable-static
make -j $(nproc) install

cd /home
git clone https://github.com/coin-or/Ipopt.git
cd /home/Ipopt
git checkout origin/stable/3.14

./configure --prefix=/home/build --with-hsl-cflags="$(pkg-config --cflags coinhsl)" --with-hsl-lflags="$(pkg-config --libs coinhsl)" --with-lapack="$(pkg-config --libs lapack)" --disable-linear-solver-loader --disable-sipopt --disable-java --enable-static
make -j $(nproc) install

cp --recursive /home/build /mydrive/build