!** Used  modules
module tranalmod
!*  Atomic parameters    
     real*8, allocatable :: &
          SX(:),SY(:),SZ(:), &   ! atom coordinates
          OX(:),OY(:),OZ(:), &   ! temporary coordinates
          WX(:),WY(:),WZ(:), &   ! local coordinates relative COM
          VX(:),VY(:),VZ(:), &   ! velocities
          Q(:)                   ! charges
!* Molecular parameters
     real*8, allocatable :: &
          X(:),Y(:),Z(:),   &    ! molecular COM coordinates
          PX(:),PY(:),PZ(:), &   ! molecular momenta
          QX(:),QY(:),QZ(:)      ! molecular angular momenta
!* Properties
     real*8, allocatable :: &
          CHARGE(:),MASS(:),SUMMAS(:)
     character(4), allocatable :: NM(:)
     character(6), allocatable :: NAME(:)
!* Reference arrays
     integer, allocatable :: &
          ITYPE(:),NNUM(:),NSITE(:),NUMR(:),ITS(:), &
          IADDR(:),ISADR(:),ISADDR(:),ITM(:),LIST(:)
     logical, allocatable :: LMMOL(:)
   end module tranalmod
   
   module tcfmod
     real*4, allocatable :: CC1(:),CC2(:),CC3(:),CC4(:),CC5(:),CC6(:)
     real*8, allocatable :: &
          DPX(:),DPY(:),DPZ(:),UPX(:),UPY(:),UPZ(:), &
          RMX(:,:,:), &
          CFA(:),CFV(:),CP1(:),CP2(:),CU1(:),CU2(:), &
          CVX(:),CVY(:),CVZ(:),CAX(:),CAY(:),CAZ(:),CCF(:)
     integer, allocatable :: N1(:),N2(:),INX(:)
     logical LMOL1,LMOL2
     character*6 NAMECF(12)
     character*2 NMM*2
     integer NCALLS,NSTEG,NSTCF,NOM,ITCF(18)
   end module tcfmod
