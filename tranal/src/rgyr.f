*   Part of MDynaMix package v.5.3 
*
*   Block TRANAL_BASE - set up and read MD trajectories
C     Computation of gyration radius (Rg) and end-to-end distance (Ree) of a molecule
C
C     Input parameters:
C     RMAX : Maximum distant to compute distribution of Rg and Ree
C     NA   : Number of bins for histograms of Rg and Ree
C     ITYP : Type of the molecules to compute Rg and Ree
C     IS1, IS2 : sites determining two ends of the molecule      
C     FILOUT: output file      
      program RGYR
      include "tranal.h"
      integer*8, allocatable :: IR2(:),IRG(:)
      character*32 FILOUT
      namelist /RGLIST/ NA,RMAX,ITYP,IS1,IS2,FILOUT
      call SETTRAJ
      read(*,RGLIST)
      allocate(IR2(NA),IRG(NA))
      IEND=0
      do I=1,NA
         IR2(I)=0
         IRG(I)=0
      end do
      SR2=0.
      SRG=0.
      SR22=0.
      SRG2=0.
      do while(IEND.eq.0)
        call READCONF(IEND)
        call GETCOM
        if(IPRINT.ge.6)write(*,'(a,f12.3)')' time',FULTIM
C     cycle over molecules
        do IMOL=IADDR(ITYP)+1,IADDR(ITYP+1)
           ISHIFT=ISADDR(ITYP)+(IMOL-1)*NSITS(ITYP)
           ISIT1=ISHIFT+IS1
           ISIT2=ISHIFT+IS2
           DX=SX(ISIT1)-SX(ISIT2)
           DY=SY(ISIT1)-SY(ISIT2)
           DZ=SZ(ISIT1)-SZ(ISIT2)
*           call PBC(DX,DY,DZ)
           R2=sqrt(DX**2+DY**2+DZ**2)
           if(R2.le.RMAX)then
             NR2=NA*R2/RMAX+1
             IR2(NR2)=IR2(NR2)+1
           end if 
           SR2=SR2+R2
           SR22=SR22+R2**2
           RGS=0.
           do JS=1,NSITS(ITYP)
             JST=JS+ISADR(ITYP) 
             JSIT=ISHIFT+JS
             DX=SX(JSIT)-X(IMOL)
             DY=SY(JSIT)-Y(IMOL)
             DZ=SZ(JSIT)-Z(IMOL)
             call PBC(DX,DY,DZ)
             RGS=RGS+MASS(JST)*(DX**2+DY**2+DZ**2)
           end do
           RG=sqrt(RGS/SUMMAS(ITYP))
           if(IPRINT.ge.7)write(*,*)IMOL,R2,RG
           if(RG.le.RMAX)then
             NRG=NA*RG/RMAX+1
             IRG(NRG)=IRG(NRG)+1
           end if  
           SRG=SRG+RG
           SRG2=SRG2+RG**2
         end do  
      end do
      R2=SR2/(IAN*NSPEC(ITYP))
      RG=SRG/(IAN*NSPEC(ITYP))
      R22=SR22/(IAN*NSPEC(ITYP))
      RG2=SRG2/(IAN*NSPEC(ITYP))
      DR2=sqrt(R22-R2**2)
      DRG=sqrt(RG2-RG**2)
      write(*,*)IAN,' configurations analysed'
      open(unit=4,file=FILOUT,status='unknown')
      write(4,'(a)')'# Molecule ',NAMOL(ITYP)
      write(4,'(2(a,f12.3))')'# End-to-end distance:',R2,'  +/- ',DR2
      write(4,'(2(a,f12.3))')'# Radius of gyration: ',RG,'  +/- ',DRG
      write(4,'(a)')'# Distributions:'
      write(4,'(a)')'#   R         <R2>         Rg  '
      do I=1,NA
         RR=(I-0.5)*RMAX/NA
         FACT=NA*1./(IAN*NSPEC(ITYP))
         write(4,'(3f12.3)')RR,FACT*IR2(I),FACT*IRG(I)
      end do   
      stop
      end
      
      
