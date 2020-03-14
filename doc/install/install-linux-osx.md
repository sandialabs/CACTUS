Compiling CACTUS for Linux/Mac
================================================

1. Uncompress the directory and put in desired location. From a terminal `cd`
   into `CACTUS/make`.
2. Create the executable with Make (from a terminal):  

   For pgf95 compiler:  
   ```
   make
   ```  

   For gfortran compiler:  
   ```
   make -f Makefile.gfortran
   ```  

   For ifort compiler:  
   ```
   make -f Makefile.ifort
   ```  

   The executable will be called `cactus` and be located in the `bin` directory.
   Alternative makefiles must be used to enable OpenMP acceleration. For
   gfortran with OpenMP:
   ```
   make -f Makefile.gfortran.omp
   ```  

   OpenMP makefiles for other compilers have not yet been written, but it should
   be easy to create them by adding the appropriate compiler flags.  

3. The compilation can be tested with the bundled regression tests.
   Navigate into the `test/RegTest/`  and run
   ```
   PATH=$PATH:../../bin/ pytest ./runreg.py
   ```
   This requires a Python installation with the `pytest` package installed.
