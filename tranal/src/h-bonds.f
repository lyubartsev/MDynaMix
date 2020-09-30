C    TRANAL v 5.3
C    Hydrogen bonds statistics
C    
C    Parameters: 
C     FOUT - name of the output file (char) 
C     RMAX - maximum donor-acceptor distance  (default 3.4 Å)
C     HANG - maximum HB angle (degrees, default 30)
C     NHB - number of different types of HBs
C     ACC - list of acceptors (as site nummer, integer array of length NHB)
C     DON - list of donors (as site nummer, integer array of length NHB)
C     HYD - list of hydrogens (as site nummer, integer array of length NHB)
      program HBONDS
      parameter(NHBMAX=1000)
      include "tranal.h"
      integer*8 ACC(NHBMAX),DON(NHBMAX),HYD(NHBMAX),ICNT(NHBMAX)
      character*64 FOUT
      namelist /HBOND/FOUT,RMAX,HANG,NHB,ACC,DON,HYD
      call SETTRAJ
      RMAX=3.4
      HANG=30.
      read(*,HBOND)
      ICONF=0
      ICNT=0
      RR2M=RMAX**2
      HANGR=HANG*3.14159256/180.
      open(unit=4,file=FOUT,status='unknown')
      do while(IEND.eq.0)
        call READCONF(IEND)
* cycle over HB types         
        do I=1,NHB
           ISA=ACC(I)   ! sites
           ISD=DON(I)
           ISH=HYD(I)
           ITA=ITS(ISA)      ! types
           ITD=ITS(ISD)
           if(ITS(ISH).ne.ITD)then
              write(*,*)' !!! Error: !!!'
              write(*,*)' HB type:                   ',I
              write(*,*)' Molecule type of donor:    ',ITD
              write(*,*)' Molecule type of hydrogen: ',ITH
              write(*,*)' Must be the same !!! '
              stop
           end if
*     cycles over molecules
           do MolD=1,NSPEC(ITD)
              IAD=ISADDR(ITD)+(MolD-1)*NSITS(ITD)+ISD-ISADR(ITD) 
              IAH=ISADDR(ITD)+(MolD-1)*NSITS(ITD)+ISH-ISADR(ITD) 
              do MolA=1,NSPEC(ITA)
*  computation of DA distance
                IAA=ISADDR(ITA)+(MolA-1)*NSITS(ITA)+ISA-ISADR(ITA)
                DX=SX(IAA)-SX(IAD)
                DY=SY(IAA)-SY(IAD)
                DZ=SZ(IAA)-SZ(IAD)
                RR2=DX**2+DY**2+DZ**2
**                write(*,'(3i7,2x,3a4)')IAD,IAH,IAA,
**     + NM(NSITE(IAD)),NM(NSITE(IAH)),NM(NSITE(IAA))
                if(RR2.le.RR2M)then
* Computation of angle
                  HX=SX(IAH)-SX(IAD)
                  HY=SY(IAH)-SY(IAD)
                  HZ=SZ(IAH)-SZ(IAD)
                  RH2=HX**2+HY**2+HZ**2
                  SCP=HX*DX+HY*DY+HZ*DZ
                  ANG=acos(SCP/sqrt(RH2*RR2))
***                  HBK=RMAX-HANG*ANG**2
***                  RR=sqrt(RR2)
***                  if(RR.le.HBK)then
                  if(ANG < HANGR)then
!     this is HB!
                    ICNT(I)=ICNT(I)+1
                  end if   
                end if     
              end do
           end do   
        end do   
      end do
      write(4,*)' Statistics of hydrogen bonds:'
      do I=1,NHB
         ANUM=ICNT(I)*1./IAN
         AANUM=ANUM/NSPEC(ITS(ACC(I)))
         ADNUM=ANUM/NSPEC(ITS(DON(I)))
         write(4,'(i4,a,a4,a,a4,a,a4,f12.3,a,2(f9.5,a))')
     + I,':  ',NM(DON(I)),'-',NM(HYD(I)),'-> ',NM(ACC(I)),
     + ANUM, ' H-bonds, ', AANUM, ' per acceptor', ADNUM, ' per donor'
      end do   
      stop
      end
