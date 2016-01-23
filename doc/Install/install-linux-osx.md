Installation of CACTUS for Linux or Mac machines
================================================

1. Uncompress the directory and put in desired location. From a terminal `cd` into the correct directory.
2. Create the executable with Make (from a terminal):  

   For pgf95 compiler:  
   ```
   make -f make/Makefile
   ```  

   For gfortran compiler:  
   ```
   make -f make/Makefile.gfortran
   ```  

   For ifort compiler:  
   ```
   make -f make/Makefile.ifort
   ```  

   The executable will be called `cactus` and be located in the `bin` directory.
   Alternative makefiles must be used to enable OpenMP acceleration. For
   gfortran with OpenMP:  
   ```
   make -f make/Makefile.gfortran.omp
   ```  
   
   OpenMP makefiles for other compilers have not yet been written, but it should
   be easy to create them by adding the appropriate compiler flags.  

3. Add the executable's location to your path so that you can call it by simply
   typing `cactus` on the terminal instead of the complete path. In Linux this
   is generally done in the `~/.bashrc` file.  

4. Move into the `test/RegTest` directory and run the regression tests by
   executing `runreg.py` with the path to the CACTUS executable as an argument
   (or just the name of the executable if you added it to your path):  
   ```
   python ./runreg.py ../../bin/cactus
   ```
   or
   ```
   python ./runreg.py cactus
   ```  

   The output should say `No differences` for all three tests. NOTE: You will
   need Python installed to run the regression tests.
