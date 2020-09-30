C    Calculation of average distance and variance between atom pairs
      program DIST
      include "tranal.h"
      parameter (MAXPAR=1000)
      real*8 RRS(MAXPAR),RRS2(MAXPAR)
      integer IS(MAXPAR),JS(MAXPAR)
      character*32 FOUT
      call SETTRAJ
      STR=TAKESTR(5,IE)
      read(STR,*)FOUT
      open(unit=8,file=FOUT,status='unknown')
      STR=TAKESTR(5,IE)
      read(STR,*)NPAR
      if(NPAR.gt.MAXPAR)NPAR=MAXPAR
      do I=1,NPAR
         STR=TAKESTR(5,IE)
         read(STR,*)IS(I),JS(I)
         RRS(I)=0.
         RRS2(I)=0.
      end do
      IEND=0
      do while(IEND.eq.0)
        call READCONF(IEND)
        if(IEND.ne.0)go to 10
        do I=1,NPAR
          II=IS(I)
          JJ=JS(I)
          IIT=ITS(II)
          JJT=ITS(JJ)
          if(IIT.ne.JJT)then
             write(*,*)' pair ',I,' is of different type'
             stop
          end if
          do IM=1,NSPEC(IIT)
            IT=ISADDR(IIT)+(IM-1)*NSITS(IIT)+IS(I)-ISADR(IIT)
            JT=ISADDR(IIT)+(IM-1)*NSITS(IIT)+JS(I)-ISADR(IIT)
            DX = SX(IT)-SX(JT)
            DY = SY(IT)-SY(JT)
            DZ = SZ(IT)-SZ(JT)
            call PBC(DX,DY,DZ)
            RR2 = DX**2+DY**2+DZ**2
            RR=sqrt(RR2)
            RRS(I)=RRS(I)+RR
            RRS2(I)=RRS2(I)+RR2
          end do
        end do  ! (NPAR
        write(*,'(a,f12.3,3f11.3)')
     +  ' time',FULTIM,RR,RR2
      end do
 10   write(*,*)IAN,' configurations analysed'
      write(8,*)IAN,' configurations analysed'
      write(8,*)'  Average distances and variations:'
      write(8,*)'      pair     <R>       <dR>        '
      IMLK=IAN*NSPEC(IIT)
      do I=1,NPAR
        IMLK=NSPEC(ITS(IS(I)))*IAN
        RR=RRS(I)/IMLK
        DR=RRS2(I)/IMLK-RR**2
        write(8,'(2i6,2f12.4)')IS(I),JS(I),RR,DR
      end do
      stop
      end
