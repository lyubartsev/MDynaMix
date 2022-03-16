!     
! File:   xdrwrap.f90
! Author: amirzoev
!
! Created on March 21, 2016, 3:29 PM
!

module xdrwrap
#ifdef XDR    
    use xdr    
    
    implicit none
    
    type(xtcfile) :: xtc, xtc_out
    type(trrfile) :: trr
#endif
    contains
#ifdef XDR
    subroutine xtc_open(fname, NumAtoms,ierr, mode)
        character(*), intent(in):: fname
        integer(4),intent(out):: NumAtoms, ierr
        character(1), optional :: mode
        
        if (present(mode) .and. (mode .eq. 'w')) then ! writing file
            call xtc_out % init(fname, 'w')
            ierr = xtc_out % stat        
        else ! reading file
            call xtc % init(fname)
            NumAtoms = xtc % NATOMS
            ierr = xtc % stat        
        endif                
    end subroutine

    subroutine trr_open(fname,NumAtoms,ierr)
        character(*), intent(in):: fname
        integer(4),intent(out):: NumAtoms, ierr
        
        call trr % init(fname)
        NumAtoms = trr % NATOMS
        ierr = trr % stat                
    end subroutine
    
    subroutine xtc_read_tranal(BOXL,BOYL,BOZL,NOP,SX,SY,SZ,FULLTIME,ierr)
        integer(4),intent(out):: ierr
        real(8), intent(out):: BOXL,BOYL,BOZL, SX(:), SY(:), SZ(:), FULLTIME
        integer(4) N,NOP
        
        N=xtc%NAtoms
        if(N.ne.NOP)then
           write(*,*)'!!! Wrong number of atoms:'
           write(*,*)'    In trajectory:              ',N
           write(*,*)'    According to the input file:',NOP
           stop
        end if
        call xtc % read
        ierr = xtc % stat
        if (ierr.eq.0) then
            SX(1:N)=xtc%pos(1,1:N)*10.d0
            SY(1:N)=xtc%pos(2,1:N)*10.d0
            SZ(1:N)=xtc%pos(3,1:N)*10.d0
            BOXL=xtc%box(1,1)*10.d0
            BOYL=xtc%box(2,2)*10.d0
            BOZL=xtc%box(3,3)*10.d0
            FULLTIME=xtc%time
!            write(*,*) 'Precision:', xtc % prec
        endif        
    end subroutine
    
    subroutine xtc_read(BOX, Ri, FULLTIME,ierr)
        integer(4),intent(out):: ierr
        real(4), intent(out):: BOX(3),  FULLTIME
        real(8), intent(out) :: Ri(:,:)
        integer(4) N
        N=xtc%NAtoms
        call xtc % read
        ierr = xtc % stat
        if (ierr.eq.0) then
            Ri(1:3, 1:N)=xtc%pos(1:3,1:N)*10.d0
            BOX=[xtc%box(1,1), xtc%box(2,2), xtc%box(3,3)]*10.d0
            FULLTIME=xtc%time
!            write(*,*) 'Precision:', xtc % prec
        endif        
    end subroutine

    
    subroutine trr_read(BOXL,BOYL,BOZL,NOP,SX,SY,SZ,VX,VY,VZ,FULLTIME,ierr)
        integer(4),intent(out):: ierr
        real(8), intent(out):: BOXL,BOYL,BOZL, SX(:), SY(:), SZ(:),VX(:), VY(:), VZ(:), FULLTIME
        integer(4) N,NOP
        
        N=trr%NAtoms
        if(N.ne.NOP)then
           write(*,*)'!!! Wrong number of atoms:'
           write(*,*)'    In trajectory:              ',N
           write(*,*)'    According to the input file:',NOP
           stop
        end if
        
        call trr % read
        ierr = trr % stat
        if (ierr.eq.0) then
            SX(1:N)=trr%pos(1,1:N)*10.d0
            SY(1:N)=trr%pos(2,1:N)*10.d0
            SZ(1:N)=trr%pos(3,1:N)*10.d0
            VX(1:N)=trr%vel(1,1:N)*10.d0
            VY(1:N)=trr%vel(2,1:N)*10.d0
            VZ(1:N)=trr%vel(3,1:N)*10.d0

            BOXL=trr%box(1,1)*10.d0
            BOYL=trr%box(2,2)*10.d0
            BOZL=trr%box(3,3)*10.d0
            FULLTIME=trr%time
        else
            call xtc % close
        endif
    end subroutine

    subroutine xtc_write(natoms, BOX, R_vec, FULLTIME, ierr)
        integer(4),intent(out):: ierr
        real(4), intent(in):: R_vec(:,:)
        real(4), intent(in):: FULLTIME, box(3)
        integer(4), intent(in) :: NAtoms

        real(4) dbox(3,3)
        dbox = reshape([box(1), 0.0, 0.0, 0.0, box(2), 0.0, 0.0, 0.0, box(3)],[3,3])
        call xtc_out % write(natoms, int(fulltime/1000), fulltime, dbox , R_vec, 1000.0)
        
        ierr = xtc_out % stat
    end subroutine
    
    subroutine xtc_close(mode)
        character(1), optional :: mode
        if (present(mode) .and. (mode .eq. 'w')) then
            call xtc_out % close
        else
            call xtc % close
        endif
    end subroutine
    
#endif    
end module xdrwrap
