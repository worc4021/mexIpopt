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
CC=icx-cl CFLAGS="-nologo -Qstd=c89"  ./configure --prefix=[your/binary/prefix] --enable-msvc --enable-static --disable-shared
make -j $(nproc)
make install
```

1. Now, we build the [coinhsl](https://github.com/coin-or-tools/ThirdParty-HSL)/[mumps](https://github.com/coin-or-tools/ThirdParty-Mumps) and which ever else solvers you might have access to. I did not encounter issues with this part, hence a variation of 
```
CC=icx-cl FC=ifx F77=ifx ./configure --prefix=[your/binary/prefix] --enable-msvc --enable-static --disable-shared --with-metis-cflags="-I[your/binary/prefix]/include/coin-or/metis" --with-metis-lflags="[your/binary/prefix]/lib/coinmetis.lib"
make -j $(nproc)
make install
```
sets you up for both of them. 
> Notice, that the installed targets for `pkg-config` don't work satisfactorily (`pkg-config --static ipopt` e.g. returns nothing for me), hence we have to help it out a little. The `CFLAGS` for the solvers are straight forward copy-past of `pkg-config --cflags coinmetis` (or coinmumps, coinhsl etc.). The `LFLAGS` are slightly more manual: Again, call `pkg-config --libs coinmetis` (or whichever) and replace the libraries it wants to link with the absolute path of the `.lib`. 

1. Next up, we build the core of ipopt. Very similar stuff, we deactivate all the dynamic loading parts and components we don't need/implement.

```
FC=ifx CC=icx-cl CXX=icx-cl F77=ifx ./configure --prefix=[you/binary/prefix] --enable-static --disable-shared --disable-linear-solver-loader --enable-msvc --disable-java --disable-sipopt --with-hsl-cflags="-I[you/binary/prefix]/include/coin-or/hsl -I[you/binary/prefix]/include/coin-or/metis" --with-hsl-lflags="[you/binary/prefix]/lib/libcoinhsl.lib [you/binary/prefix]/lib/libcoinmetis.lib /c/PROGRA~2/Intel/oneAPI/mkl/latest/lib/mkl_intel_lp64.lib /c/PROGRA~2/Intel/oneAPI/mkl/latest/lib/mkl_sequential.lib /c/PROGRA~2/Intel/oneAPI/mkl/latest/lib/mkl_core.lib" --with-mumps-cflags="-I[you/binary/prefix]/include/coin-or/mumps -I[you/binary/prefix]/include/coin-or/metis" --with-mumps-lflags="[you/binary/prefix]/lib/libcoinmumps.lib [you/binary/prefix]/lib/libcoinmetis.lib /c/PROGRA~2/Intel/oneAPI/mkl/latest/lib/mkl_intel_lp64.lib /c/PROGRA~2/Intel/oneAPI/mkl/latest/lib/mkl_sequential.lib /c/PROGRA~2/Intel/oneAPI/mkl/latest/lib/mkl_core.lib"
make -j $(nproc)
make install
```

The logic is the same, we copy paste the `.lib` files to where the `-lfoo` points and pass the `CFLAGS` through.

> Regarding the tests that ipopt ships with: I found that they did not work on mingw despite working fine when executed in powershell directly. You can of course build them `make test` and then try them from windows, but I'd expect them to fail when run by `make` itself.

# Building the mex

This part is actually rather straight forward if you got this far. Make sure the submodules are all present. If you have a similar setup with MKL etc you should be able to simply run

```
cmake --preset oneApi
cmake --build --preset oneApi-build
ctest --preset oneApi-test
```

As mentioned above, once you have a mex file, it depends on nothing outside your matlab installation. I.e. you can copy it onto any other Matlab installation of compatible Matlab release (i.e. using a later version is fine).

## Notes:
- When rebuilding, I found that often I needed to clear the previous build folder if any configuration of upstream libraries changed.