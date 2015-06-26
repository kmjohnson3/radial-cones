This repository contains source code to generate the radial-cones trajectory described in MRM [ publication link to be updated]. 

There are two sets of code.

# Generation Code #
This matlab based code will generate the base cone trajectory.  See Generation/radial_cones_example.m

# Optimization code #
This C++ based code optimizing the radial trajectory utilizing and outpute file "RadialCones.struct" from the matlab code. This code exports a set of rotation matrix of 3 x 3 x #shots. This code requires Armadillo (Matrix library), FFTW (dft), Blitz++ (array library), HDF5_cpp (i/o).  This code was run in a Ubunutu 14.04 base machine with the libraries configures as below:

```
#!c++


"acml-5-3-1-gfortran-64bit.tar" (support for armadillo)
./install-acml-5-3-1-gfortran-64bit.sh -accept -installdir=${HOME}/local/acml/  

"hdf5-1.8.14.tar"
./configure --prefix=${HOME}/local/ --enable-fortran --enable-cxx --enable-static --disable-shared  

"armadillo-4.600.2.tar"
./configure -DCMAKE_INSTALL_PREFIX=${HOME}/local/   

"fftw-3.3.4.tar"
./configure --enable-float --enable-threads --enable-openmp --prefix=${HOME}/local/

"blitz-0.10.tar"
./configure --prefix=${HOME}/local/



```