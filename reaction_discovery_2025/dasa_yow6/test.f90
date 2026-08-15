program debug_file_read

    implicit none
    integer :: numbr, irr, i, ios
    real :: vbarTMP
    integer :: qq2(14)
    character(len=100) :: foundp
    logical :: eof_reached

    ! Assuming 'foundp' is a filename, replace it with the actual filename if needed
    foundp = 'foundp.txt'
    
    open(14, file=trim(foundp), status='old', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error opening file', ios
        stop
    end if

    eof_reached = .false.
    do i = 1, 20  ! Loop for more iterations to test beyond the problematic point
        write(*,*) 'beginning of do, iteration:', i
        
        ! Check for end-of-file before reading
        inquire(unit=14, end=eof_reached)
        if (eof_reached) then
            write(*,*) 'EOF reached, exiting loop'
            exit
        end if

        ! Read numbr and vbarTMP
        write(*,*) 'before read'
        read(14,*,iostat=irr) numbr, vbarTMP
        if (irr /= 0) then
            write(*,*) 'Error reading numbr and vbarTMP, iostat:', irr
            close(14)
            exit
        end if
        write(*,*) 'after read, iostat = ', irr
        write(*,*) "foundp", numbr
        
        ! Read first set of 5 integers
        read(14,*,iostat=irr) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
        if (irr /= 0) then
            write(*,*) 'Error reading first set of 5 integers, iostat:', irr
            close(14)
            exit
        end if
        
        ! Read second set of 5 integers
        read(14,*,iostat=irr) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
        if (irr /= 0) then
            write(*,*) 'Error reading second set of 5 integers, iostat:', irr
            close(14)
            exit
        end if
        
        ! Read last set of 4 integers
        read(14,*,iostat=irr) qq2(11), qq2(12), qq2(13), qq2(14)
        if (irr /= 0) then
            write(*,*) 'Error reading last set of 4 integers, iostat:', irr
            close(14)
            exit
        end if
        
        write(*,*) 'read complete, iteration:', i
    end do

    close(14)
    write(*,*) 'Program completed successfully'

end program debug_file_read

