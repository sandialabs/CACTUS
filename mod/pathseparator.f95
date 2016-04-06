module pathseparator

    implicit none

    character(1) :: path_separator

contains
    
    subroutine get_path_separator()
        ! get_path_separator() : Get the path delimiter ('/' for Linux, '\' for Windows)
        character(len=99999) :: path

        call get_environment_variable('PATH',path)

        path_separator = path(1:1)

    end subroutine get_path_separator

end module pathseparator