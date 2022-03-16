C    Tranal - package for rajectory analysis, MDynaMix v.> 5.0
C    Compute Ramachandra plot for 2 torsion angles
C    
      program TORSION2
      include "tranal.h"
      integer, allocatable :: IT(:),JT(:),KT(:),LT(:)
      integer, allocatable :: NPOP(:,:)
      character FTOR*64,OTYP*16
      call SETTRAJ
*  Input of data - after "TRAJ" section, 
*  positional input with optional commentaries
*  
*  1) File name of the output file and output type ("2D" or "gnu")
      IE=0
      STR=TAKESTR(5,IE)
      read(STR,*,end=12,err=12)FTOR,OTYP
      if(OTYP(1:1).eq.'t'.or.OTYP(1:1).eq.'T')then
        IOTYP=2
      else
        IOTYP=1
      end if
      go to 13
 12   read(STR,*)FTOR
      IOTYP=1
 13   continue
*  2) number of different pairs of torsions and number of bins in [-180:180]
      STR=TAKESTR(5,IE)
      read(STR,*)MANG,NAA
      allocate(NPOP(NAA,NAA))
      NADT=MANG*2               ! total
      allocate(IT(NADT),JT(NADT),KT(NADT),LT(NADT))
*   3) List of torsions - 4 site numbers for each torsion
*      All "phi" torsions, then all "psi" torsions
      do I=1,NADT                   
	STR=TAKESTR(5,IE)
  	READ(STR,*)IT(I),JT(I),KT(I),LT(I)    ! which
      end do
      if(IPRINT.ge.6)then 
	write(*,*)' Compute distributions over torsion angles:'
	do I=1,MANG
          J=2*I-1
          JTYP=ITS(J)
	  write(*,*)NM(NSITE(IT(J))),NM(NSITE(JT(J))),
     + NM(NSITE(KT(J))),NM(NSITE(LT(J))),'mol typ',JTYP
          J=J+1
          JTYP=ITS(J)
	  write(*,*)NM(NSITE(IT(J))),NM(NSITE(JT(J))),
     + NM(NSITE(KT(J))),NM(NSITE(LT(J))),'mol typ',JTYP
          write(*,*)'--------------------------------------'
	end do               
      end if
      do I=1,NAA
	do J=1,NAA
	  NPOP(I,J)=0 
	end do
      end do
      MTYP=ITS(IT(1))
      do I=1,NADT
        II=ITS(IT(I))
        JJ=ITS(JT(I))
        KK=ITS(KT(I))
        LL=ITS(LT(I))
        if(II.ne.MTYP.or.JJ.ne.MTYP.or.KK.ne.MTYP.or.LL.ne.MTYP)then
           write(*,*)'Different molecular types in the list of torsion.'
           write(*,*)'Check the input'
           stop
        end if
      end do
      write(*,*)' Molecule type is ',MTYP
*  start analysis
      IEND=0
      do while(IEND.eq.0)
        if(IEND.ne.0)go to 101
        call READCONF(IEND)
        if(IPRINT.ge.6)write(*,*)'Proceeding time ',FULTIM
	  do IMOL=1,NSPEC(MTYP)         !  over molecules of this type
	  do N=1,MANG
            do IN=1,2
            ITOR=2*(N-1)+IN
            I         = IT(ITOR)  
            J         = JT(ITOR)
            K         = KT(ITOR)
            L         = LT(ITOR)
	    IST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+I
	    JST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+J
	    KST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+K
	    LST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+L
        ax        = sx(jst)-sx(ist)
        ay        = sy(jst)-sy(ist)
        az        = sz(jst)-sz(ist)
        bx        = sx(kst)-sx(jst)
        by        = sy(kst)-sy(jst)
        bz        = sz(kst)-sz(jst)
        cx        = sx(lst)-sx(kst)
        cy        = sy(lst)-sy(kst)
        cz        = sz(lst)-sz(kst)
*
        call PBC(AX,AY,AZ)
        call PBC(BX,BY,BZ)
        call PBC(CX,CY,CZ)
*
        ab        = ax*bx + ay*by + az*bz
        bc        = bx*cx + by*cy + bz*cz
        ac        = ax*cx + ay*cy + az*cz
        at        = ax*ax + ay*ay + az*az
        bt        = bx*bx + by*by + bz*bz
        ct        = cx*cx + cy*cy + cz*cz

        axb       = (at*bt)-(ab*ab)
        bxc       = (bt*ct)-(bc*bc)
*
        fnum      = (ab*bc)-(ac*bt)
        den       = axb*bxc
        if(den.gt.0.0d0) then 
*
          den       = dsqrt(den)         
Cosine of angle:
          co        = fnum/den
          CO        = DMIN1(CO, 1.D0)
          CO        = DMAX1(CO,-1.D0)
* Sign of angle:
      signum    = ax*(by*cz-cy*bz)+ay*(bz*cx-cz*bx)+az*(bx*cy-cx*by)
* Value of angle:
          arg       = dsign(dacos(co),signum)
	  else
	    arg=0.d0
	  end if    
          if(IN.eq.1)then
	    IPHI=0.5*(arg+PI)*NAA/PI+1
          else
	    IPSI=0.5*(arg+PI)*NAA/PI+1
          end if
	  if(IPRINT.ge.8)write(*,'(7i5,f10.3)')
     +N,IN,IMOL,IST,JST,KST,LST,ARG*180/PI
	  end do
          NPOP(IPHI,IPSI)=NPOP(IPHI,IPSI)+1
	  end do
	end do
      end do
 101  write(*,*)IAN,' configurations analysed'
      open(unit=17,file=ftor,status='unknown')
      FAC=NAA**2*1./(IAN*MANG*NSPEC(MTYP))
      if(IOTYP.eq.1)then
         write(17,'(a)')'# distribution over torsion angles'
         do I=1,NAA
	   PHI=-180+(I-0.5)*360/NAA 
           do J=1,NAA
	     PSI=-180+(J-0.5)*360/NAA 
             write(17,'(2f10.4,f12.5)')PHI,PSI,NPOP(I,J)*FAC
           end do
         end do
      else
         do I=1,NAA
           write(17,'(1000f12.5)')(NPOP(I,J)*FAC,J=1,NAA)
         end do
      end if
      close(17)
      write(*,*)' Results are written in file ',FTOR
      stop
      end
