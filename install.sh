#!/bin/bash
apt-get update
apt-get -y install liblapack-dev
apt-get -y install mpich
apt-get -y install gfortran
apt-get -y install libopenmpi-dev
apt-get -y install openmpi-bin
apt-get -y install gcc-multilib
apt-get -y install libgsl0ldbl
apt-get -y install libgsl-dev
apt-get -y install libcgal-dev
apt-get -y install cmake
apt-get -y install vim
apt-get -y install git
apt-get -y install unp
# install ONLY cuda toolkit according to http://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#abstract (use download .deb variant)
# download anaconda 64 bit linux
# download google chrome 64 bit
# download Eigen3
