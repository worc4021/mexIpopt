# Introduction 

On Linux and on Macos getting an ipopt mexfile is rather straight forward - on Windows a completely different story. Hence, most of this description is to get you going on Windows.

## Prerequisites 

On Windows, I tested this with the oneAPI compilers through [Mingw64](https://www.mingw-w64.org/). I use `zsh` on mingw, if all this means nothing to you and you have mingw and the icx compiler installed simply make a batch script with 

```
CALL "C:\\Program Files (x86)\\Intel\\oneAPI\\setvars.bat" intel64 vs2022
CALL "C:\\msys64\\msys2_shell.cmd" -defterm -here -use-full-path -no-start -mingw64 -shell zsh
```

and adapt your paths/shell if necessary. You should also have a _current_ version of [cmake](https://cmake.org/) installed. You will need 3.21 or later.

# Preparing the additional requirements

We need some linear solvers before we can build ipopt itself. Nowadays, the MKL's pardiso implementation is also supported by ipopt but you might want to others of their supported linear solvers.

Due to the nature of this being a project for a mex implementation of ipopt, I will only show how to use a statically linked build. This means that once you have a mex, all you need to copy into your matlab project is the mex file, no additional library load paths or shared objects to copy.

1. First off, we build [metis](https://github.com/coin-or-tools/ThirdParty-Metis). After you run `get.Metis` you will want to auto configure the build. Metis does not use headers to declare functions, and your compiler will throw errors about undeclared symbols. To circumvent this you need to instruct your compiler to run in c89 mode. I.e. your configuration looks something like this:
```
CC=icx-cl CFLAGS="-nologo -Qstd=c89" ADD_CFLAGS=-MT  ./configure --prefix=[your/binary/prefix] --enable-msvc --enable-static --disable-shared
make -j $(nproc)
make install
```
We will use the prefix folder to do most of the heavy lifting for us. Namely the `pkg-config` `.pc` files `autoconfigure` installs for us. Subsequent build steps can read them to configure themselves.
1. Now, we build the [coinhsl](https://github.com/coin-or-tools/ThirdParty-HSL)/[mumps](https://github.com/coin-or-tools/ThirdParty-Mumps) and which ever else solvers you might have access to. I did not encounter issues with this part, hence a variation of 
```
CC=icx-cl FC=ifx F77=ifx ADD_CFLAGS=-MT ADD_FCFLAGS=-MT ADD_FFLAGS=-MT PKG_CONFIG_PATH=[prefix]/lib/pkgconfig ./configure --prefix=[prefix] --enable-msvc --enable-static --disable-shared
make -j $(nproc)
make install
```
sets you up for both of them. 
> Notice that since we are only building static objects we do not need the linker information from `pkg-config`, we continue to use the windows frontend of icx.

1. Next up, we build the main ipopt. Very similar stuff, we deactivate all the dynamic loading parts and components we don't need/implement, however, to build the tests as well as to detect the available componentns we will need to link against the libraries we have just built. Hence, we require that `-L/c/myprefix` etc makes sense to our compiler, this is not the case for the msvc frontend `icx-cl` hence we switch to `icx` for this.

```
FC=ifx CC=icx CXX=icx F77=ifx ADD_CFLAGS=-MT ADD_FCFLAGS=-MT ADD_FFLAGS=-MT ADD_CXXFLAGS=-MT PKG_CONFIG_PATH=[prefix]/lib/pkgconfig ./configure --prefix=[you/binary/prefix] --enable-msvc --enable-static --disable-shared --disable-linear-solver-loader --disable-java --disable-sipopt 
make -j $(nproc)
make install
```

> The tests will require that you also have the intel runtimes installed, otherwise the build will work but the test might fail.

Note, that if you want to build the debug symbols version you need to use `-MTd` as well as `--enable-debug` for all the configurations above.

# Building the mex

This part is actually rather straight forward if you got this far. Make sure the submodules are all present. If you have a similar setup with MKL etc you should be able to simply run

```
cmake --preset windows-intel-release-config
cmake --build --preset windows-intel-release-build
ctest --preset windows-intel-release-test
```

As mentioned above, once you have a mex file, it depends on nothing outside your matlab installation. I.e. you can copy it onto any other Matlab installation of compatible Matlab release (i.e. using a later version is fine).

## Notes:
- When rebuilding, I found that often I needed to clear the previous build folder if any configuration of upstream libraries changed.


## Building with oneApi on linux

It turned out that compiling on linux with the oneApi had a few minor challenges, I used the following:

```
CC=icx CFLAGS="-std=c90 -fpie" ./configure --prefix=[prefix] --enable-static --disable-shared
make -j $(nproc)
make install
```

again Metis was having issues with the language standard and its undefined symbols. The position independent executable was something needed for the ipopt tests themselves, i.e. for executables.

Subsequently, on linux the `pkg-config` can do most of the heavy lifting, so Mumps and the HSL compile ok with 

```
CC=icx FC=ifx F77=ifx PKG_CONFIG_PATH=[prefix]/lib/pkgconfig ADD_CFLAGS="-fpie" ADD_FCFLAGS="-fpie" ADD_FFLAGS="-fpie" ./configure --prefix=[prefix] --enable-static --disable-shared
make -j $(nproc)
make install
```

And ultimately the ipopt library itself:

```
FC=ifx CC=icx CXX=icpx F77=ifx ADD_CFLAGS="-fpie" ADD_CXXFLAGS="-fpie" ADD_FCFLAGS="-fpie" ADD_FFLAGS="-fpie" PKG_CONFIG_PATH=[prefix]/lib/pkgconfig ./configure --prefix=[prefix] --enable-static --disable-shared --disable-linear-solver-loader --disable-java --disable-sipopt
make -j $(nproc)
make test
make install
```