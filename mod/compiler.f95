#define xstr(s) str(s)
#define str(s) #s

module compiler

    ! compiler module

    implicit none

    ! Include mkl.fi include file if using Intel Fortran compiler
#if defined(__INTEL_COMPILER)
    include "mkl.fi"
#endif

    character(80) :: COMPILER_STRING, COMPILER_VER

contains
    

    subroutine print_compiler_info()

    ! print_compiler_info() : print compiler info to stdout

#if defined(__GFORTRAN__)
        COMPILER_VER = __VERSION__
        write(*,'(A)') "Compiled with GNU Fortran ", COMPILER_VER
#elif defined(__INTEL_COMPILER)
        COMPILER_VER = xstr(__INTEL_COMPILER)
        write(*,'(A,A,".",A)') "Compiled with Intel Fortran version: ", COMPILER_VER(1:2), COMPILER_VER(3:4)
#else
        write(*,'(A)') "Compiled with unknown compiler."
#endif

    ! Alert if MKL is enabled
#if defined(__INTEL_COMPILER)
!dec$if defined(__INTEL_MKL__)
        write(*,'(A)') "Intel MKL enabled."
!dec$endif
#endif
    end subroutine print_compiler_info


end module compiler
