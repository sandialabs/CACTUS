# Compiling CACTUS

## General instructions
Compilation requires CMake, which can be installed by your OS package manager, or manually (see https://cmake.org/install/)

1. Clone the repo, or download all source files to a folder.

2. In the root of the repo, create a `build/` directory, and run `cmake`/`make`
```
mkdir -p build
cd build
cmake ../ 
make
```
Support for OpenMP may be disabled by adding the `DOPENMP=OFF` flag to the `cmake` call, as below.
```
cmake -DOPENMP=OFF ../ 
```

3. The compilation can be tested with the bundled regression tests.
   Navigate into the `test/RegTest/`  and run
   ```
   PATH=$PATH:../../bin/ pytest ./runreg.py
   ```
   This requires a Python 3 installation with the `pytest` package installed.

## Windows Tips
For Windows, it is recommended to compile in MingW under MSYS2. 

1. Install MSYS2 (https://www.msys2.org/)
2. Open MingW terminal (likely located at `C:\tools\msys64\mingw64.exe`)
3. Install the following packages using `pacman`
```
pacman -Sy git make cmake mingw-w64-x86_84-gcc-fortran mingw-w64-x86_64-lapack
```
4. Follow the general install instructions above.

For MSYS builds, the executables will be statically linked for portability.
