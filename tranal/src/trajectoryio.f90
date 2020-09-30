MODULE TrajectoryIO
    implicit none   
    private
    public:: SaveFrame, OpenTraj, CloseTraj, ReadFrame
contains
    subroutine SaveFrameXMOL(stdout_, NAtom, BOX, Names, R_vec, time)
        implicit none
        integer(4), intent(in)   :: stdout_,  NAtom
        character(*), intent(in) :: Names(:)
        real(4), intent(in)      :: BOX(3), time
        real(8), intent(in)      :: R_vec(:,:)

        integer(4) iAtom

        !  First Line - number of atoms
        write(stdout_,*) NAtom
    !  Second commentary line of XMOL format contains time in fs or current MC-step
        write(stdout_,'(a7,f14.0,a,3f10.3)')' after ',time,' BOX:', box        ! for MC
        do iAtom = 1, NAtom
            write(stdout_,'(a,3f14.5)') Names(iAtom), R_vec(1:3, iAtom)
        end do               
    end subroutine

#ifdef XDR
    subroutine SaveFrameXTC(stdout_, NAtom, BOX, R_vec, time)
        use xdrwrap
        implicit none
        integer(4), intent(in)   :: stdout_,  NAtom
        real(4), intent(in)      :: BOX(3), time
        real(8), intent(in)      :: R_vec(:,:)

        integer(4) iAtom, ierr
        
        ierr = 0
        call xtc_write(natom, BOX*0.1, real(R_vec*0.1, 4), TIME*0.001, ierr)
        if (ierr.ne.0) stop "Error while writing output trajectory file"
    end subroutine
#endif    

    subroutine SaveFrame(fmt, stdout_, NAtom, BOX, Names, R_vec, time)            
        use strings, only: str_cln
        implicit none
        integer(4), intent(in)   :: stdout_,  NAtom
        character(*), intent(in) :: Names(:)
        real(4), intent(in)      :: BOX(3), time
        real(8), intent(in)      :: R_vec(:,:)
        character(*), intent(in) :: fmt
        
        select case (str_cln(fmt))
            case ('XMOL')            
                call SaveFrameXMOL(stdout_, NAtom, BOX, Names, R_vec, time)
#ifdef XDR            
            case ('XTC')        
                call SaveFrameXTC(stdout_, NAtom, BOX, R_vec, time)                
#endif
            case default
                stop "Unknown file format: "
        end select
    end subroutine
    
    subroutine OpenTraj(fmt, filename, stdout_, mode, ierr)
        use strings, only: str_cln
        use xdrwrap
        implicit none
        integer(4), intent(in)   :: stdout_
        character(*), intent(in) :: fmt, filename
        character(1), intent(in) :: mode
        integer(4), intent(out) :: ierr
        integer(4) dummy
        
        ierr = 0
        
        select case (str_cln(fmt))
#ifdef XDR              
            case('XTC')                                  
                select case (str_cln(mode))
                    case ('R','W')
                        call xtc_open(filename, dummy, ierr, mode)
                    case default
                        stop "Unsupported file reading mode. Must be either 'R' or 'W' "                   
                end select                    
#endif                            
            case ('XMOL')
                select case (str_cln(mode))
                    case ('R')
                        open(unit=stdout_, file=filename, status='OLD', iostat=ierr, action='READ') 
                    case ('W')
                        open(unit=stdout_, file=filename, IOSTAT=ierr)
                    case default
                        stop "Unsupported file reading mode. Must be either 'R' or 'W' "                   
                end select    

            case default
                write(*,*) 'Error OpenTraj - unsupported trajectory format:', str_cln(fmt)
                stop 657
        end select
        if (ierr.ne.0) stop "Error while opening output trajectory file"
    end subroutine
    
    subroutine CloseTraj(fmt, stdout_, mode)
    use strings, only: str_cln
    use xdrwrap
    implicit none
    character(*), intent(in):: fmt
    character(1), intent(in):: mode
    integer(4), intent(in) :: stdout_

        select case(str_cln(fmt))
#ifdef XDR                        
            case ('XTC')
                call xtc_close(mode)
#endif                            
            case ('XMOL')
                close(stdout_)
            case default
                stop "Unsupported file format: "                              
        end select        
    end subroutine
    
    subroutine ReadFrame(fmt, stdin_, NAtom, BOX, Names, R_vec, time, ierr)
        use precision
        use strings, only: str_cln
        implicit none
        integer(4), intent(in)   :: stdin_
        integer(4), intent(out) :: NAtom
        integer(4), intent(out):: ierr        
        character(*), intent(in):: fmt
        character(atom_name_len), intent(inout) :: Names(:)
        real(4), intent(out)      :: BOX(3), time
        real(8), intent(out)      :: R_vec(:,:)
        
        select case(str_cln(fmt))
#ifdef XDR                        
            case ('XTC')
                call ReadFrameXTC(BOX, R_vec, time, ierr)
#endif                            
            case ('XMOL')
                call ReadFrameXMOL(stdin_, NAtom, BOX, Names, R_vec, time, ierr)
            case default
                stop "Unsupported file format: "                              
        end select   
    end subroutine
    
    subroutine ReadFrameXMOL(stdin_, NAtom, BOX, Names, R_vec, time, ierr)
        use precision
        implicit none
        integer(4), intent(in)   :: stdin_
        integer(4), intent(out) :: NAtom
        integer(4), intent(out):: ierr        
        character(atom_name_len), intent(inout) :: Names(:)
        real(4), intent(inout)      :: BOX(3), time
        real(8), intent(out)      :: R_vec(:,:)
        character(input_text) line
        integer(4) iAtom,pos_box, pos_time1, pos_time2 !, NAtom_

        !  First Line - number of atoms
        read(stdin_,*, iostat=ierr) NAtom
        if (is_iostat_end(ierr)) then
                write(*,*) 'End of the  file  reached.'
                return
        else if (ierr.ne.0) then
            write(*,*)'ERROR: while reading from file '
            stop 10
        endif        
        read(stdin_, '(a)', iostat=ierr) line
!         write(*,*) line
!        if (NAtom.ne.0) then ! i.e. if we knew number of atoms
!            if (NAtom_.ne.NAtom) then
!                write(*,*) 'Error ReadFrameXMOL: Incorrect number of atoms in the frame:', NAtom_, ' Shall be:', NAtom
!                stop 659
!            endif
!        endif
!         write(*,*) line
        pos_box = index(line, 'BOX:')
        if (pos_box.ne.0) then
            read(line(pos_box+4:), *) BOX(:)        
        else
            BOX = 0.0
        endif
    !  Second commentary line of XMOL format may contain time in fs or current MC-step, let's try to detect it        
        pos_time1 = index(line, 'after')
        if (pos_time1.ne.0) pos_time1 = pos_time1 + 5
        pos_time2 = index(line, 'fs')
        if ((pos_time1.eq.0) .and. (pos_time2.eq.0)) then   ! didn't detect from first try, make another attempt assuming LAMMPS.xyz format        
                pos_time1 = index(line, 'Timestep:')
                if (pos_time1.ne.0) pos_time1 = pos_time1 + 9
                pos_time2 = len_trim(line)
        endif
!        write(*,*) pos_time1, pos_time2, line(pos_time1+5:pos_time2)
        if ((pos_time1.ne.0) .and. (pos_time2.ne.0)) then           
            read(line(pos_time1:pos_time2), *) time        
        endif
       
        do iAtom = 1, NAtom
            read(stdin_,*,iostat=ierr) Names(iAtom), R_vec(1:3, iAtom)
        end do                       
    end subroutine
#ifdef XDR
    subroutine ReadFrameXTC(BOX, R_vec, time, ierr)       
        use xdrwrap
        use precision
        implicit none
!        integer(4), intent(in)   :: stdin_
!        integer(4), intent(inout) :: NAtom
        integer(4), intent(out):: ierr        
!        character(atom_name_len), intent(inout) :: Names(:)
        real(4), intent(inout) :: BOX(3), time        
        real(8), intent(out)      :: R_vec(:,:)
        call xtc_read(BOX, R_vec, time, ierr)
        if (ierr.ne.0) then
            write(*,*) 'Error ReadFrameXTC: Error while reading trajectory.', ierr
            stop 656
        endif
    end subroutine
#endif    
    
END MODULE TrajectoryIO