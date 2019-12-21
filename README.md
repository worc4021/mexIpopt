# Introduction 
1. Build coin-hsl library
1. Build the mkl blas and lapack libraries.
1. Build ipopt as a library
1. Generate mex interface.
1. Done

# Build COIN-HSL solver
Run 
```
./configure --prefix=${binFolder}
```
in the coinhsl repo to set up all relevant make files etc. 
Then run 
```
make -j8
make install
```
Now there should be a coinhsl library and archive in the prefix directory.

# Build BLAS and LAPACK
Run 
```
make libuni export=blas_example_list name=libblas
make libuni export=lapack_example_list name=liblapack
mv *.dylib ${binFolder}
```
This creates MKL exports for us to link against.


# Build Ipopt
Get the release you want from the Ipopt repo. Then run
```
./configure --prefix=${binFolder} --with-hsl-cflags="-I${binFolder}/include" --with-hsl-lflags="-L${binFolder}/lib -lcoinhsl" --with-lapack="-L${binFolder}/lib -llapack" CC="/usr/bin/xcrun gcc" CXX="/usr/bin/xcrun g++" F77="/usr/bin/xcrun gfortran" --disable-linear-solver-loader
```
The linear solver loader for some reason does only work when we have access to all HSL solvers, since we don't we turn it off.

Presumably, there is some left over garbage in the compiler paths pointing to no longer existing fortran runtimes. I was not able to fix this. So I ended up simply modifying the make files. Soz.

The run 
```
make -j8
make install
```
done.

# Build mex interface
Run the make file in the repository.
