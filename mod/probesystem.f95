module probesystem

    ! probe module for computing and writing velocities at a specified location

    use configr
    use pathseparator

    implicit none

    ! probe output options
    integer :: ProbeFlag                     ! flag for enabling probes (1 for probes, 0 for none)
    integer :: ProbeOutIntervalTimesteps     ! interval between probe writes
    integer :: ProbeOutStartTimestep         ! timestep at which to start writing probe data
    integer :: ProbeOutEndTimestep           ! timestep at which to end writing probe data

    ! probe output header flag
    logical :: ProbeHeadersWritten = .False.

    type ProbeType

        real          :: p(3)      ! probe location
        real          :: vel(3)    ! probe velocity (at a single instant)
        real          :: vel_fs(3) ! probe free-stream velocity (at a single instant)
        character(80) :: output_filename

    end type ProbeType

    type(ProbeType), allocatable :: probes(:)         ! array to hold probes
    integer                      :: nprobes           ! number of probes
    integer, parameter           :: probe_iounit = 28 ! IO unit for probes


contains

    subroutine compute_probe_vel(probe)

        ! compute_probe_vel() : computes the induced velocity and free-stream velocity
        !   at the probe location and stores it to the probe variables

        type(ProbeType), intent(inout) :: probe
        integer                        :: ygcErr ! error flag for CalcFreeStream

        ! Calculate wall and wake induced velocities at grid locations
        Call CalcIndVel(NT,ntTerm,NBE,NB,NE, &
            probe%p(1),probe%p(2),probe%p(3), &
            probe%vel(1),probe%vel(2),probe%vel(3))

        ! Calculate free stream velocities at grid locations
        Call CalcFreestream(probe%p(1),probe%p(2),probe%p(3), &
            probe%vel_fs(1),probe%vel_fs(2),probe%vel_fs(3), &
            ygcErr)

    end subroutine compute_probe_vel


    subroutine read_probes(probespec_filename)
        ! read_probes() : reads in the coordinates of probes from a file
        !   Coordinates should be formatted as
        !   x1 y1 z1
        !   x2 y2 z2
        !   ...
        !   EOF

        integer i
        integer, parameter           :: iounit = 20
        character(len=*), intent(in) :: probespec_filename ! Probe specification filename

        ! open file for reading
        open(unit=iounit, form='formatted', file=probespec_filename)

        ! read number of probes
        read(iounit,*) nprobes

        ! allocate for probes
        allocate(probes(nprobes))

        ! read probe coordinates
        do i=1,nprobes
            read(iounit,*) probes(i)%p(1), probes(i)%p(2), probes(i)%p(3)
        end do

        ! close file
        close(unit=iounit)

    end subroutine read_probes


    subroutine write_probe_headers(probe_output_path)

        ! write_probe_headers() : writes the headers for the probe output files to the specified
        !   probe_output_path

        integer                     :: probe_num
        character(10000), parameter :: probe_file_header='X/R (-),Y/R (-),Z/R (-)'
        character(10000), parameter :: probe_data_header='Normalized Time (-),U/Uinf (-),V/Uinf (-),W/Uinf (-),Ufs/Uinf (-),Vfs/Uinf (-),Wfs/Uinf (-)'
        character(80)               :: probe_num_str
        character(80)               :: probe_output_path

        ! write header files
        do probe_num=1,nprobes

            ! assemble the probe filename in the pattern
            !   [run_name]_probe_[probe_num].csv
            write(probe_num_str,'(I5.5)') probe_num
            probes(probe_num)%output_filename = 'probe_'//trim(probe_num_str)//'.csv'

            ! open/write header/close file
            open(probe_iounit, file=trim(adjustl(probe_output_path))//path_separator//probes(probe_num)%output_filename)

            ! write probe location
            write(probe_iounit,'(A)') trim(probe_file_header)
            write(probe_iounit,'(E13.7,",",$)') probes(probe_num)%p(1)
            write(probe_iounit,'(E13.7,",",$)') probes(probe_num)%p(2)          ! probe location
            write(probe_iounit,'(E13.7,","  )') probes(probe_num)%p(3)

            ! write probe data header
            write(probe_iounit,'(A)') trim(probe_data_header)

            ! close file
            close(probe_iounit)

        end do

    end subroutine write_probe_headers


    subroutine write_probes(probe_output_path)

        ! write_probes() : Computes the velocity for all probes and writes the data to the
        !    corresponding probe output file.

        type(ProbeType) :: probe
        integer         :: probe_num
        character(80)   :: probe_output_path

        ! Write the headers to the probe output files if they haven't been written already
        if (ProbeHeadersWritten .eqv. .False.) then
            call write_probe_headers(probe_output_path)
            ProbeHeadersWritten = .True.
        end if

        ! Write the probe data
        do probe_num=1,nprobes

            ! get the current probe
            probe = probes(probe_num)

            ! open file for writing
            open(probe_iounit, file=trim(adjustl(probe_output_path))//path_separator//probe%output_filename, POSITION='append')

            ! compute velocity at probe location
            call compute_probe_vel(probe)

            ! write the position, velocity, and free-stream velocity
            write(probe_iounit,'(E13.7,",",$)') TimeN               ! Normalized simulation time (t*Uinf/Rmax)
            write(probe_iounit,'(E13.7,",",$)') probe%vel(1)
            write(probe_iounit,'(E13.7,",",$)') probe%vel(2)        ! induced velocities
            write(probe_iounit,'(E13.7,",",$)') probe%vel(3)
            write(probe_iounit,'(E13.7,",",$)') probe%vel_fs(1)
            write(probe_iounit,'(E13.7,",",$)') probe%vel_fs(2)    ! freestream velocities
            write(probe_iounit,'(E13.7)'      ) probe%vel_fs(3)    ! Dont suppress carriage return on last column

            ! close file
            close(probe_iounit)

        end do

    end subroutine write_probes

end module probesystem