# CACTUS
CACTUS (**C**ode for **A**xial and **C**ross-flow **TU**rbine **S**imulations), developed at Sandia National Laboratories, is a turbine simulation code based on a free wake vortex method. This repository is a fork of the original CACTUS repository which includes a number of additional features and improvements, including:
- Output of wake data (Cartesian and filament data)
- Modeling of general wall geometries
- OpenMP acceleration of induced velocity calculations

Contributors:
Phillip Chiu (pchiu@sandia.gov)

### Installation & Compilation
Installation and compilation instructions can be found in `Doc/Install`

### Directory Structure
- `Airfoil_Section_Data`: Contains airfoil data tables.
- `bin`: Executables are placed here.
- `CreateGeom`: Contains MATLAB geometry creation scripts
- `DAKOTA`: Contains DAKOTA Drivers and examples
- `Doc`: Documentation: User's Manual, Install instructions, DAKOTA-CACTUS Manual, Papers
- `Makefiles`: Contains MakeFiles for different compilers
- `mod`: Source Code -- modules & utilities
- `src`: Source Code
- `Test`: Contains several test cases

### Post-processing
Tools for post-processing data from CACTUS simulations are available in the [CACTUS-tools ](https://github.com/whophil/CACTUS-tools]) repository.

### References
For details about the development of CACTUS, please see
- Murray, J., and Barone, M., “The Development of CACTUS, a Wind and Marine Turbine Performance Simulation Code,” _49th AIAA Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition_, Reston, Virginia: American Institute of Aeronautics and Astronautics, 2011, pp. 1–21.
