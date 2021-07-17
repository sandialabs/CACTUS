module pathseparator

    implicit none

    character(1) :: path_separator

contains
    
    subroutine get_path_separator()

#if defined(PATHSEP)
        path_separator = PATHSEP
#else
        ! fallback
        path_separator = '/'  
#endif
    end subroutine get_path_separator

end module pathseparator