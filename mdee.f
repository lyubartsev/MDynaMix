*
*   Molecular Dynamics Expanded Ensemble
*
*   1. Start matters
*
	subroutine MDEE_INIT
	include "common_inc.f"
	if(TASKID.eq.MAST)then
	  write(*,*)' Molecular Dynamics - Expanded Ensemble'
	  if(LL2MC)then
	    IDEE=-1 
            write(*,*)' L2MC algorithm for changes of subensembles'
	  end if  
	  write(*,*)' Molecules to annihilate:'
	  end if
	NDEL=0	  
	if(.not.LRST)then
          IWLI = 1             ! num of WL iteration
	  IPASS = 0            ! current num of EE steps in the current sweep
	  MPASS = 0            ! num of EE steps when pass condition was reached
	end if
	do I=1,NTYPES
	  if(IEE(I).eq.0)then
	    NDEL=NDEL+NSPEC(I)
	    if(TASKID.eq.MAST)write(*,'(4x,a32,i6,a)')
     +      NAMOL(I),NSPEC(I),' mol'
	  end if
	end do
	if(NDEL.eq.0)then
	  write(*,*)
     + '!!! MDEE was specified but no EE molecules were defined'
	  write(*,*)'    going back to normal MD'
	  write(*,*)
	  LMDEE=.false.
	  return
	end if
        do IE=1,NE
          if(IPRINT.ge.7)write(*,*)EC(IE),EE(IE)
C  balancing factors EE are given per molecule in the input file
	  EE(IE)=EE(IE)*NDEL
        end do
C  starting subensemble    
	if(.not.LRST)ME=ME0
C  setup interactions
	call RECINT
        call DIFF_LRCOR
C  setting zero arrays
	if(.not.LRST)then
	  do I=1,NENS
	    do J=1,NENS
	      IWALK(I,J)=0
	    end do
	    do J=1,LHIST
	      ITE(I,J)=0
	    end do
	  end do
	end if	
	if(LAUTOFACT)then
          WL_CORR=WL_INI*WL_FAC**(IWLI-1)*NDEL
	  if(IWLI.le.NWLI)then
	    if(TASKID.eq.MAST)then
	      write(*,*)
     +' Wang-Landau algorithm for optimization of balancing factors'
	      write(*,*)' WL iteration no ',IWLI
	      write(*,*)' WL factor       ',WL_CORR
	   end if
	  else
	    if(TASKID.eq.MAST)write(*,*)
     +' Continue EE simulation with fixed balancing factors'
	    WL_CORR=0.d0
	  end if
	else
	  WL_CORR=0.d0
	end if
	NNH=NE+3
	return
	end
*	
*======================= RECINT ===============================
*
*   2. Changing interactions parameters in EE
*
*   Rule: LJ_eps is scaled by alpha**RW1
*         Charges scaled by alpha**PW1
*
      subroutine RECINT
      include "common_inc.f"
      SPEI=0.
C  molecules of first type
      do ITYP=1,NTYPES
C  whether scale this molecule in EE
	if(IEE(ITYP).eq.0)then
	  ESCI=EC(ME)**PW1 
	else
	  ESCI=1.d0
	end if   
C   cycle over sites (I)
        DO I       =  ISADR(ITYP)+1,ISADR(ITYP+1)
	  IS=INBT(I)
C   molecules of the second type
	  do JTYP=1,NTYPES
C   total LJ EE scaling factor (if I or J is EE molecule)
	    ESCL=ESCI
	    if(IEE(JTYP).eq.0)ESCL=EC(ME)**PW1
C  cycle over sites (J)
            DO J       =  ISADR(JTYP)+1,ISADR(JTYP+1)
	      JS=INBT(J)
C  I,J: interaction sites
C  IS,JS - interaction types
C  The rule (subject to optional change)
*  Must be consistent with similar fragment in DIFEN and DIF_LRCOR below
	      A6LJ(IS,JS)=ESCL*A6LJ0(IS,JS)
	      B12LJ(IS,JS)=ESCL*B12LJ0(IS,JS)
            END DO! OF J
          END DO
	end do
*  Scale charge
C   cycle over atoms (ISP)
        DO ISP     =  ISADDR(ITYP)+1,ISADDR(ITYP+1)
	  I = NSITE(ISP)
	  Q(ISP)=ESCI*CHARGE(I)
	end do
*  recomputing self-energy term in Ewald / Fourier
	SPET=0.
	do I=1,NSITS(ITYP)
	  IS=ISADR(ITYP)+I
  	  QII=COULF*(ESCI*CHARGE(IS))**2
          SPET         = SPET-NSPEC(ITYP)*ALPHAD*QII/dsqrt(PI)
	end do
	SPET=SPET/NUMTASK
	SPEEI(ITYP)=SPET
	SPEI=SPEI+SPET
      end do
*  shift of distance
      RAEX = BEE*(1.-EC(ME)**PW2)
      return
      end
*
*   3. Change of ensemble
*======================= CHENS ================================
*
	subroutine CHENS
	include "common_inc.f" 
	real*8 BUFF(4),BUF2(4)
C  Decide up or down
	if(TASKID.eq.MAST)then
	  if(LL2MC)then 
	    ME1=ME+IDEE
	  else
	    if(rand().lt.0.5)then
	      ME1=ME-1
	    else
	      ME1=ME+1
	    end if
	  end if   
	end if  
C  report to other nodes
	call BICAST(ME1,MAST)
C  Increase counter of ensemble change attempts
	time20=cputime(0.d0)   
	ICHE=ICHE+1
	IPASS=IPASS+1
	if(ME1.lt.1.or.ME1.gt.NE)then
	   if(LL2MC)then
	     ITRAN=0
	     IDEE=-IDEE 
	   end if
	   go to 90
	end if   
C   Compute change of energy due to reciprocal part of Ewald  
	call DIFFUR(ME1,EFN,EFO)
C   Compute change of energy due to other interactions
	call DIFEN(ME1,EEN,EEO)
C  Sum up the total energy and send to "master" node
	BUFF(1)=EFN
	BUFF(2)=EFO
	BUFF(3)=EEN
	BUFF(4)=EEO
	call SUMMA(BUFF,BUF2,4,MAST)
	EFN=BUF2(1)
	EFO=BUF2(2)
	EEN=BUF2(3)
	EEO=BUF2(4)
        if(IPRINT.ge.7)write(*,'(a,i3,4f15.4)')
     +  ' node ',TASKID,EEN*BETA,EEO*BETA,EFN*BETA,EFO*BETA
	if(TASKID.eq.MAST)then
	  DEF=EFN-EFO
	  DEE=EEN-EEO  
	  DE=(DEE+DEF+EELRC(ME1)-EELRC(ME))*BETA-EE(ME1)+EE(ME)
	  ENNN=(EEN+EFN+EELRC(ME1))*BETA
	  ENNO=(EEO+EFO+EELRC(ME))*BETA
	  if(IPRINT.ge.6)write(*,'(a5,I5,a4,I5,a4,f12.4,3x,3f11.3)')
     + 'from ',ME,' to ',ME1,' DE=',DE,ENNN,ENNO,EE(ME1)-EE(ME)
C   MC transition - Metropolis rule
	  if(exp(-DE).lt.rand())then
	    ITRAN=0
	    IDEE=-IDEE 
	  else
	    ITRAN=1
	  end if
	end if
 	call BICAST(ITRAN,MAST)
C  Successful transition
	if(ITRAN.eq.1)then
	  ME=ME1
	  ICHU=ICHU+1
   	  call RECINT
C   Redefine interactions
C  Add to transition matrix
          do I=1,ME-1
            if(IWALK(I,ME).le.IWALK(ME,I))
     +      IWALK(I,ME)=IWALK(I,ME)+1
          end do
          do I=ME+1,NE
            if(IWALK(I,ME).lt.IWALK(ME,I))
     +      IWALK(I,ME)=IWALK(I,ME)+1
          end do 
C  Checking pass condition
	  if(MPASS.eq.0.and.(ME.eq.1.or.ME.eq.NE))then
	    if(IWALK(1,NE).ge.NPASS.and.IWALK(NE,1).ge.NPASS)
     +        MPASS=IPASS*(1.+0.5/NPASS)
	  end if
	end if
C  Decrease of the balancing factors in case of WL
 90	if(LAUTOFACT.and.IWLI.le.NWLI)then
          EE(ME)=EE(ME)-WL_CORR
C  Completing WL iteraction
	  if(MPASS.gt.1.and.IPASS.gt.MPASS)then
	    do I=1,NE
              EE(I)=EE(I)-EE(NE)
              do J=1,NE
                IWALK(I,J)=0
              end do
            end do
            WL_CORR=WL_INI*WL_FAC**IWLI*NDEL
	    HIST(3,NHIST)=1.
	    if(IPRINT.ge.4)then
	      write(*,'(I4,a,I12)')
     + IWLI,' WL iteraction completed. Num. steps:',MPASS*NCEE
	      write(*,*)' New balancing factors:'
	      write(*,'(10(10f8.3))')(EE(IC),IC=1,NE)
	    end if
	    if(IWLI.lt.NWLI)then
	      if(IPRINT.ge.4)write(*,*)' New WL increment: ',WL_CORR
	    else
	      if(IPRINT.ge.4)write(*,*)
     + ' Continue with fixed balancing factors'
	      WL_CORR=0.
	    end if		
	    IHIST=NHIST+1
	    IWLI = IWLI+1
	    IPASS=0
	    MPASS=0
	  end if
        endif
 100	time2=time2+cputime(time20)
	return
	end                     
*
*===================== DIFFUR ===========================================
*
C  4. DIFFUR
C  Compute difference of energy due to reciprocal part of the Ewald sum
C
      subroutine DIFFUR(ME1,EFN,EFO)
      include "common_inc.f"
      DATA TWOPI/6.2831853072D0/
*
      timef0	= cputime(0.d0)
*
      ALP2 = ALPHAD**2
      CFUR=4.d0*PI*COULF/VOL     
*
      EFN     = 0.D0
      EFO	= 0.d0 
*
      ALEEO = EC(ME)**PW1
      ALEEN = EC(ME1)**PW1
*   4.1 Reciprocal space Evald
      IBEG = NUMTASK-TASKID 
      do IKV=IBEG,NKV,NUMTASK
 	RXV=RKX(IKV)/BOXL
	RYV=RKY(IKV)/BOYL
	RZV=RKZ(IKV)/BOZL
        RK2=RXV**2+RYV**2+RZV**2
        SCSO=0.
        SSNO=0.
	SCSN=0.
	SSNN=0.
        do I=1,NSTOT
	  IC=NSITE(I)
	  ITYP=ITYPE(I)
	  SCP=RXV*SX(I)+RYV*SY(I)+RZV*SZ(I)
          SSC=sin(SCP)
	  CSC=cos(SCP)
	  QO = Q(I)
	  if(IEE(ITYP).eq.0)then
	    QN=CHARGE(IC)*ALEEN
	  else
	    QN=Q(I)
	  end if
          SSNO=SSNO+QO*SSC
          SCSO=SCSO+QO*CSC
          SSNN=SSNN+QN*SSC
          SCSN=SCSN+QN*CSC
        end do
	AKC = CFUR*exp(-0.25d0*RK2/ALP2)/RK2
        EFO=EFO+(SSNO**2+SCSO**2)*AKC
        EFN=EFN+(SSNN**2+SCSN**2)*AKC
      end do
*
*	FACTR=  0.001*ENERF/FNOP
*	if(IPRINT.ge.7)write(*,'(a,4f12.5)')' FURtot:',FACTR*EFN,FACTR*EFO
*   4.2. Self-interaction Ewald
*   4.2.1.  Constant term
CC	write(*,*)' 1 ',EFO*BETA,EFN*BETA
      ABC = ALPHAD/DSQRT(PI)
      do ITYP=1,NTYPES
	if(IEE(ITYP).eq.0)then
	  ESCN = ALEEN
	  ESCO = ALEEO
	  IBEG = NUMTASK - TASKID
          do II=IBEG,NSITS(ITYP),NUMTASK
	    IS=ISADR(ITYP)+II
	    EFO = EFO - COULF*(ESCO*CHARGE(IS))**2*ABC*NSPEC(ITYP)
	    EFN = EFN - COULF*(ESCN*CHARGE(IS))**2*ABC*NSPEC(ITYP)
	  end do
	end if
      end do
CC	write(*,*)' 2 ',EFO*BETA,EFN*BETA
*  4.2.2  Intramolecular term  (see ETERMS in forces.f)
C  must substract contribution of reciprocal Ewald for intramolecular 
C  electrostatic interactions
      do ITYP=1,NTYPES
*  only for EE molecules
	if(IEE(ITYP).eq.0.and.NSITS(ITYP).gt.1)then
*  over pairs of sites of this type
          ISBEG       = ISADR  (ITYP)+1
          ISEND       = ISADR  (ITYP+1)
          do IS=ISBEG+TASKID,ISEND-1,NUMTASK
	    do JS=IS+1,ISEND
*  over molecules
	      do IM = 1,NSPEC(ITYP)
		ISP = ISADDR(ITYP)+(IM-1)*NSITS(ITYP)+IS-ISBEG+1 
		JSP = ISADDR(ITYP)+(IM-1)*NSITS(ITYP)+JS-ISBEG+1 
   	        QIJO = COULF*Q(ISP)*Q(JSP)
	        QIJN = COULF*ALEEN**2*CHARGE(IS)*CHARGE(JS)
                DX          = SX(ISP)-SX(JSP)
                DY          = SY(ISP)-SY(JSP)
                DZ          = SZ(ISP)-SZ(JSP)
	        call PBC(DX,DY,DZ)
                R2          = DX**2+DY**2+DZ**2
                R2I         = 1.D0/R2
                R1          = DSQRT(R2)
                R1I         = R1*R2I
                ALPHAR      = ALPHAD*R1
*  	          call RERFC(ALPHAR,ERFC,EXP2A)
                TTT         = 1.D0/(1.D0+B1*ALPHAR)
                EXP2A       = DEXP(-ALPHAR**2)
                ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A
	        ERFR = (ERFC-1.d0)*R1I
	        EFO = EFO + QIJO*ERFR
                EFN = EFN + QIJN*ERFR 
	      end do
	    end do
	  end do  ! of IS
	end if  
      end do ! of ITYP
C	write(*,*)' 3 ',EFO*BETA,EFN*BETA
      timef	= timef+cputime(timef0)
      return
      end  
*
*====================== DIFEN ========================================
*           
*       5. DIFEN
*       Difference in energy due to LJ and real-space Ewald
*                                                          
	subroutine DIFEN(ME1,EEN,EEO)
	include "common_inc.f"
        timel0	= cputime(0.d0)
	EEN=0.d0
	EEO=0.d0         
	ELJO=0.
	ESTO=0.
	ELJN=0.
        ESTN=0.
*  alpha for old and new subensembles
	EEALO=EC(ME)
	EEALN=EC(ME1)
	EEALO1=EEALO**PW1
	EEALN1=EEALN**PW1
*  EE distance shift (old and new)
	RSHO = BEE*(1.-EEALO**PW2)  
	RSHN = BEE*(1.-EEALN**PW2)  
*  this is medium-range interactions: sum over all pairs
	DO INP      = 1,MBLN
	  ISP=iabs(NBL1(INP))
	  JSP=NBL2(INP)
	  ITYP=ITYPE(ISP)
	  JTYP=ITYPE(JSP)
* select relevant to EE
	  if(IEE(ITYP)*IEE(JTYP).ne.1)then
	    ISB =NSITE(ISP)
  	    JSB =NSITE(JSP)
	    IMOL=NNUM(ISP)
	    JMOL=NNUM(JSP)
*  not within the same molekule
	    if(IMOL.ne.JMOL)then
              ISN =INBT(ISB)           !   non-bonded type
              JSN =INBT(JSB)
  	      A6O =A6LJ(ISN,JSN)
	      B12O=B12LJ(ISN,JSN)
	      QIJO=Q(ISP)*Q(JSP)*COULF
	      if(IEE(ITYP).eq.0)then
	        ESCI=EEALN1 
	      else
	        ESCI=1.d0
	      end if
	      ESCL = ESCI
	      if(IEE(JTYP).eq.0)then
		ESCL=ESCI*EEALN1
	        ESCI = EEALN1
	      end if
	      A6N = ESCI*A6LJ0(ISN,JSN)
	      B12N = ESCI*B12LJ0(ISN,JSN)
              QIJN       =  ESCL*CHARGE(ISB)*CHARGE(JSB)*COULF
*
              DX           = SX(ISP)-SX(JSP)
              DY           = SY(ISP)-SY(JSP)
              DZ           = SZ(ISP)-SZ(JSP)
              call PBC(DX,DY,DZ)
*  
              RR           = DX**2+DY**2+DZ**2
	      R1  = DSQRT(RR) 
              R1O           = R1 + RSHO
              R1N           = R1 + RSHN
	      R1OI = 1./R1O
	      R1NI = 1./R1N
*  Electrostatic interaction
              ALPHAR       = ALPHAD*R1
	      TTT          = 1.D0/(1.D0+B1*ALPHAR)
              ERFCR=((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)
     &              *TTT*EXP(-ALPHAR**2)
              EESO          = QIJO*ERFCR*R1OI
              EESN          = QIJN*ERFCR*R1NI
*  LJ interaction
              R6IO          = R1OI**6
              R6IN          = R1NI**6
              ESRO          = (     A6O+      B12O*R6IO)*R6IO 
              ESRN          = (     A6N+      B12N*R6IN)*R6IN
              EEN=EEN+EESN+ESRN
	      EEO=EEO+EESO+ESRO
	      ELJN=ELJN+ESRN
	      ESTN=ESTN+EESN
	      ELJO=ELJO+ESRO
	      ESTO=ESTO+EESO
	    end if !  IMOL.ne.JMOL
C	if(IPRINT.ge.7.and.ESCL.ne.1.d0)write(*,'(2i5,2f12.5,2e14.6)')
C     +ISP,JSP,ENERF*0.001*(ESRO+EESO),ENERF*0.001*(ESRN+EESN),A6O,QIJO
*
	end if
      END DO! OF INP					
      timel	= timel+cputime(timel0) 
*    Local forces 
      times0	= cputime(0.d0)
      DO INP      = 1,MBSH
	ISP=iabs(NBS1(INP))
	JSP=NBS2(INP)
	ITYP=ITYPE(ISP)
	JTYP=ITYPE(JSP)
* select relevant to EE
	if(IEE(ITYP)*IEE(JTYP).ne.1)then
	  ISB =NSITE(ISP)
  	  JSB =NSITE(JSP)
	  IMOL=NNUM(ISP)
	  JMOL=NNUM(JSP)
*  not within the same molecule
	  if(IMOL.ne.JMOL)then
              ISN =INBT(ISB)           !   non-bonded type
              JSN =INBT(JSB)
  	      A6O =A6LJ(ISN,JSN)
	      B12O=B12LJ(ISN,JSN)
	      QIJO=Q(ISP)*Q(JSP)*COULF
	      if(IEE(ITYP).eq.0)then
	        ESCI=EEALN1 
	      else
	        ESCI=1.d0
	      end if
	      ESCL = ESCI
	      if(IEE(JTYP).eq.0)then
		ESCL=ESCI*EEALN1
	        ESCI = EEALN1
	      end if
	      A6N = ESCI*A6LJ0(ISN,JSN)
	      B12N = ESCI*B12LJ0(ISN,JSN)
              QIJN       =  ESCL*CHARGE(ISB)*CHARGE(JSB)*COULF
*
              DX           = SX(ISP)-SX(JSP)
              DY           = SY(ISP)-SY(JSP)
              DZ           = SZ(ISP)-SZ(JSP)
              call PBC(DX,DY,DZ)
*  
              RR          = DX**2+DY**2+DZ**2
	      R1  = DSQRT(RR) 
              R1O           = R1 + RSHO
              R1N           = R1 + RSHN
	      R1OI = 1./R1O
	      R1NI = 1./R1N
*  Electrostatic interaction
              ALPHAR       = ALPHAD*R1
	      TTT          = 1.D0/(1.D0+B1*ALPHAR)
              ERFCR=((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)
     &              *TTT*EXP(-ALPHAR**2)
              EESO          = QIJO*ERFCR*R1OI
              EESN          = QIJN*ERFCR*R1NI
*  LJ interaction
              R6IO          = R1OI**6
              R6IN          = R1NI**6
              ESRO          = (     A6O+      B12O*R6IO)*R6IO 
              ESRN          = (     A6N+      B12N*R6IN)*R6IN
              EEN=EEN+EESN+ESRN
	      EEO=EEO+EESO+ESRO
	      ELJN=ELJN+ESRN
	      ESTN=ESTN+EESN
	      ELJO=ELJO+ESRO
	      ESTO=ESTO+EESO
	  end if   ! of IMOL.ne.JMOL
	end if
      END DO! OF INP	     
C	write(*,'(a,4f12.3)')' DIFEN:',ELJN*BETA,ELJO*BETA,ESTN*BETA,ESTO*BETA
      times	= times+cputime(times0)
      return
      end 
*
*================ DIFF_LRCOR ======================================
*
C   Change of LJ out-cutoff corrections 
C   and self contribution from Ewald sum
C   (computed only once, or after substantial volume change)
*
      subroutine DIFF_LRCOR
      include "common_inc.f"
      RC0          = DSQRT(RSSQ)
      FAC  = 4.*PI/VOL 
*  Loop over ensembles
      do I=1,NE
        EEAL3 = EC(ME)**PW2
	EELRC(I)=0.
*  molecules of first type
        DO ITYP = 1,NTYPES
          ISBEG       = ISADR(ITYP)+1
          ISEND       = ISADR(ITYP +1)
          NSPI        = NSPEC(ITYP)
	  if(IEE(ITYP).eq.0)then
	    ESCI=EC(I)**PW1
	  else
	    ESCI=1.d0
	  end if
*  molecules of the second type
          DO JTYP = ITYP ,NTYPES
            MT          = MDX  (ITYP,JTYP)
            JSBEG       = ISADR(JTYP)+1
            JSEND       = ISADR(JTYP +1)
            NSPJ        = NSPEC(JTYP)
            FNOPIJ      = DFLOAT(NSPI*NSPJ)
            if(ITYP.eq.JTYP)FNOPIJ=0.5*FNOPIJ
            if(IEE(ITYP)*IEE(JTYP).eq.0)then
              RC = RC0 + BEE*(1.-EEAL3)
            else
              RC = RC0
            end if
	    ESCL=ESCI
	    if(IEE(JTYP).eq.0)ESCL=EC(I)**PW1
*  loop over sites in each molecule
            DO ISB = ISBEG,ISEND
              IS = INBT(ISB)           !   non-bonded type
              DO JSB = JSBEG,JSEND
		JS = INBT(JSB)
		A6 = ESCL*A6LJ0(IS,JS)
                B12=ESCL*B12LJ0(IS,JS)
		if(abs(A6).ge.1.d-50.and.abs(B12).gt.1.d-50)then
		  SI3 = sqrt(-B12/A6)
		  EP   = 0.25*A6**2/B12
		else
		  SI3=0.
		  EP=0.
		end if
                ELR  = EP*SI3*(( 1.D0/9.D0)*(SI3/RC**3)**3
     X                                - ( 1.D0/3.D0)*(SI3/RC**3))*FNOPIJ 
	        EELRC(I) = EELRC(I) + 4.*PI*ELR/VOL
	      end do
            end do
	  end do
	end do   ! do ITYP
      end do    !  over ensembles
      if(IPRINT.ge.7.and.TASKID.eq.MAST)then
	write(*,*)
     +' Out-cutoff energy corrections for different ensembles:'
	do I=1,NE
	   write(*,*)I,EELRC(I)*PERMOL
	end do
      end if
      return
      end
*
*===================== Averaging of quantities in MDEE =============
*
*
*============= MDEE_AVE  =================================================
*
      SUBROUTINE MDEE_AVE(NNN)
*
      include "common_inc.f"
	character*500 STR
	dimension PEE(NENS),FEN(NENS)
     +,TMP1(LHIST),TMP2(LHIST),TINN(NTPS*(NTPS+1)/2)
*
* !!! Attention! Any corrections should be done accordingly in 
*                parts 1 and 3 (labels 1000,3000)
*
      if(NHIST.eq.0)NHIST=1
      TOKJ           = ENERF*1.D-3
      FNOPI          = 1.D0/DFLOAT(NOP)
*
      GO TO (1000,2000,2000),NNN
*
      STOP '!!! NNN OUT OF ORDER IN GETAVR !'
*
 1000 CONTINUE 
      AVE(1,ME)=AVE(1,ME)+1
      IAV            = 2
*
*Bonds:
      AVB            = 0.D0
      M              = 0
      DO ITYP        = 1,NTYPES
      NSP            = NSPEC(ITYP)
      FNSPI          = 1.D0/DFLOAT(NSP)
      NRBS           = NRB(ITYP)
      IF(NRBS.NE.0)  THEN
      FNBI           = 1.D0/DFLOAT(NRBS)
      DO I           = 1,NRBS
      M              = M+1
      IAV            = IAV+1
      BBM            = BB(M)
      AVE(IAV,ME)        = AVE(IAV,ME)+BBM
      AWE(IAV,ME)        = AWE(IAV,ME)+BBM**2
      IAV            = IAV+1
      EBM            = EB(M)*TOKJ
      AVE(IAV,ME)        = AVE(IAV,ME)+EBM
      AWE(IAV,ME)        = AWE(IAV,ME)+EBM**2
      AVB            = AVB+EBM*FNBI
      EBOND(ITYP)    = EBOND(ITYP)+EBM
      END DO! OF I
      END IF
      END DO! OF ITYP
*Angles:
      AVA            = 0.D0
      M              = 0
      DO ITYP        = 1,NTYPES
      NSP            = NSPEC(ITYP)
      FNSPI          = 1.D0/DFLOAT(NSP)
      NRAS           = NRA(ITYP)
      IF(NRAS.NE.0)  THEN
      FNAI           = 1.D0/DFLOAT(NRAS)
      DO I           = 1,NRAS
      M              = M+1
      IAV            = IAV+1
      AAM            = AA(M)*TODGR
      AVE(IAV,ME)        = AVE(IAV,ME)+AAM
      AWE(IAV,ME)        = AWE(IAV,ME)+AAM**2
      IAV            = IAV+1
      EAM            = EA(M)*TOKJ
      AVE(IAV,ME)        = AVE(IAV,ME)+EAM
      AWE(IAV,ME)        = AWE(IAV,ME)+EAM**2
      EANGL(ITYP)    = EANGL(ITYP)+EAM
      AVA            = AVA+EAM*FNAI
      END DO! OF I
      END IF
      END DO! OF ITYP
*Torsions:
      AVT            = 0.D0
      M              = 0
      DO ITYP        = 1,NTYPES
      NSP            = NSPEC(ITYP)
      FNSPI          = 1.D0/DFLOAT(NSP)
      NRTS           = NRT(ITYP)
      IF(NRTS.NE.0)  THEN
      FNTI           = 1.D0/DFLOAT(NRTS)
      DO I           = 1,NRTS
      M              = M+1
      IAV            = IAV+1
      TTM            = TT(M)*TODGR
      AVE(IAV,ME)        = AVE(IAV,ME)+TTM
      AWE(IAV,ME)        = AWE(IAV,ME)+TTM**2
      IAV            = IAV+1
      ETM            = ET(M)*TOKJ
      AVE(IAV,ME)        = AVE(IAV,ME)+ETM
      AWE(IAV,ME)        = AWE(IAV,ME)+ETM**2
      ETORS(ITYP)    = ETORS(ITYP)+ETM
      AVT            = AVT+ETM*FNTI
      END DO! OF I
      END IF
      END DO! OF ITYP
*Other energies
	PEL=0.
      DO I           = 1,NTYPES
      NSP            = NSPEC(I)
      FNSPI          = 1.D0/DFLOAT(NSP)
      SELFPE(I)      = SELFPE(I)*TOKJ*FNSPI
      PES14 (I)      = PES14 (I)*TOKJ*FNSPI
      PSR14 (I)      = PSR14 (I)*TOKJ*FNSPI
      IAV            = IAV+1
      AVE(IAV,ME)        = AVE(IAV,ME)+PES14(I)
      AWE(IAV,ME)        = AWE(IAV,ME)+PES14(I)**2
      IAV            = IAV+1
      AVE(IAV,ME)        = AVE(IAV,ME)+PSR14(I)
      AWE(IAV,ME)        = AWE(IAV,ME)+PSR14(I)**2
	PEL	= PEL+SELFPE(I)+PES14(I)
      END DO
*
Calculate total potential energy:
*
      POTEN          = 0.D0
*Intermolecular non-bonded
      DO ITYP        =     1,NTYPES
        NSPI           = NSPEC(ITYP)
        DO JTYP        =  ITYP,NTYPES
          NSPJ           = NSPEC(JTYP)
          I              =   MDX(ITYP,JTYP)
          FNSPI = 1.D0/DFLOAT(NSPI)
          IF(NSPJ.lt.NSPI)FNSPI = 1.D0/DFLOAT(NSPJ)
          POTEN          = POTEN+PELRC(I)+POTES(I)+POTLJ(I)
          POTSEL         = PELRC(I)*TOKJ*FNSPI
          POTESE         = POTES(I)*TOKJ*FNSPI
          POTLJE         = POTLJ(I)*TOKJ*FNSPI
          PEL	= PEL+POTESE
*          IAV            = IAV+1
*         AV(IAV)        = AV(IAV)+POTSEL
*         AW(IAV)        = AW(IAV)+POTSEL**2
          IAV            = IAV+1
          AVE(IAV,ME)        = AVE(IAV,ME)+POTESE
          AWE(IAV,ME)        = AWE(IAV,ME)+POTESE**2
          IAV            = IAV+1
          AVE(IAV,ME)        = AVE(IAV,ME)+POTLJE
          AWE(IAV,ME)        = AWE(IAV,ME)+POTLJE**2
          POTES(I)       = 0.D0
          POTLJ(I)       = 0.D0
        END DO! OF JTYP
      END DO! OF ITYP
* temperature
      IAV            = IAV+1
      AVE(IAV,ME)        = AVE(IAV,ME)+TEMP
      AWE(IAV,ME)        = AWE(IAV,ME)+TEMP**2
* pot.energy
      IAV            = IAV+1
      PEX	= PE*TOKJ*FNOPI
      AVE(IAV,ME)        = AVE(IAV,ME)+PEX
      AWE(IAV,ME)        = AWE(IAV,ME)+PEX**2
*  intramolecular energy
      IAV            = IAV+1
      PINTT	= PINT*TOKJ*FNOPI
      AVE(IAV,ME)        = AVE(IAV,ME)+PINTT
      AWE(IAV,ME)        = AWE(IAV,ME)+PINTT**2
      POTEN          = POTEN*TOKJ*FNOPI+PQE*TOKJ*FNOPI
	PEL	= PEL+PQE*TOKJ*FNOPI
*  pressure
      IAV	= IAV+1
      PRES= TRYCK
      AVE(IAV,ME)	= AVE(IAV,ME)+PRES
      AWE(IAV,ME)	= AWE(IAV,ME)+PRES**2     
      IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+TRYCKM
      AWE(IAV,ME)	= AWE(IAV,ME)+TRYCKM**2
*  box sizes
      IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+BOXL
      AWE(IAV,ME)	= AWE(IAV,ME)+BOXL**2
      IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+BOYL
      AWE(IAV,ME)	= AWE(IAV,ME)+BOYL**2
      IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+BOZL
      AWE(IAV,ME)	= AWE(IAV,ME)+BOZL**2
*  translational and rotational kinetic energies
	IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+TTR
      AWE(IAV,ME)	= AWE(IAV,ME)+TTR**2
*
	IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+TROT
      AWE(IAV,ME)	= AWE(IAV,ME)+TROT**2  
*
	IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+TINT
      AWE(IAV,ME)	= AWE(IAV,ME)+TINT**2
*
	IAV	= IAV+1
      AVE(IAV,ME)	= AVE(IAV,ME)+ENNO
      AWE(IAV,ME)	= AWE(IAV,ME)+ENNO**2
*
      NNAV	= IAV
      NAVT	= NAVT+1
*
      IF(IAV.GT.NRQS) THEN
      PRINT *,'!!! INCREASE NRQS TO ',IAV
      STOP 
      END IF
*
Calculate total kinetic energy:
*
      EPOTB          = 0.D0
      EKINE          = 0.D0
      DO I           = 1,NTYPES
      EKINE          = EKINE+EKIN (I)
      EPOTB          = EPOTB+EBOND(I)+EANGL(I)+ETORS(I)
      EBOND(I)       = 0.D0
      EANGL(I)       = 0.D0
      ETORS(I)       = 0.D0
      END DO
      EKINE          = EKINE*TOKJ*FNOPI
      POTEN	= PE+PINT
*
      ETOT           = POTEN+EKINE
      ITE(ME,NHIST)=ITE(ME,NHIST)+1
      HIST(2,NHIST)=HIST(2,NHIST)+1.	
*
      RETURN
*
 2000 CONTINUE
*     
      if(NNN.eq.3.and.HIST(2,NHIST).lt.0.5d0)then
	NHIST=NHIST-1
	go to 3000
      end if
      do I=1,NE
	PEE(I)=dfloat(ITE(I,NHIST))/HIST(2,NHIST)
      end do
      if(IPRINT.ge.4)then
        write(*,'(i4,10(10f7.4/))')NHIST,(PEE(I),I=1,NE)
        write(*,'(a,3f11.3,a,f10.5)')' Box: ',BOXL,BOYL,BOZL,
     +             '   Dens: ',TOTMAS*1.d-3/(VOL*UNITL**3)
      end if
      do I=1,NE
	if(NDEL.gt.0)HIST(I+3,NHIST)=EE(I)/NDEL
      end do
*
	MNSTEP         = MSTEP       
        FNAVSI         = 1.D0/DFLOAT(NAVT)
        TTIM	= (TIM+NSTEP*DT)*1.d12
	HIST(1,NHIST)=TTIM
*     
	if(NNN.eq.3)go to 3000
        NHIST	= NHIST+1
	do I=1,NNH
	  HIST(I,NHIST)=0.d0
	end do 
	if(TASKID.eq.MAST)write(*,*)' begin ',NHIST,' iteration'
	if(NHIST.eq.IHIST)then
	  NAVT	= 0
	  if(TASKID.eq.MAST)write(*,*)' New averaging'
          TIMA	= 0.
	  call ZEROAV
          DO I        = 1,NRQS 
	    do J=1,NE
            AVE(I,J)       = 0.D0
            AWE(I,J)       = 0.D0
	    end do
          END DO! OF I
	  if(.not.LAUTOFACT)then
            do I=1,NE
              do J=1,NE
                IWALK(I,J)=0
              end do
            end do
	  end if
        end if
      if(NHIST.gt.LHIST)then
        call INSHIFT(ITE,LHIST,NE,NENS)
        call SHIFT(HIST,LHIST,NNH,NRH)
        NHIST=LHIST
      end if
      return  
 3000 continue                                                     
      if(IHIST.gt.NHIST)IHIST=NHIST 
      IBEG=1
      if(IPRINT.lt.6)IBEG=IHIST
      do I=1,NHIST
	if(HIST(3,I).eq.1.)IHIST=I
      end do
      SUMT=0.
      FACD=1.d0/sqrt(NHIST-IHIST+1.d0)
* Calculation of averages
      do I=IHIST,NHIST
        SUMT=SUMT+HIST(2,I)
      end do
*
      NST=SUMT
      if(IHIST.gt.1)then
        TM0=HIST(1,IHIST-1)
      else
        TM0=0.
	IHIST=1
      end if
	do I=1,NE
	  if(AVE(1,I).lt.1.d0)AVE(1,I)=1.d0
	  do J=1,NNAV 
	    AWE(J,I)=AVE(J,I)/AVE(1,I)
	  end do
	end do
      if(TASKID.ne.MAST)return
      PRINT "(80('-'))"
      PRINT 
     +"('#** FINAL AVERAGE QUANTITIES AFTER ',I8,' STEPS: ')",NST
      write(*,*)'  from ',TM0,'ps   to ',HIST(1,NHIST),'ps'
*     
	if(IPRINT.ge.6)then
      write(*,*)' CPU time for some procedures:'
      time=timest+cputime(time0)
      write(*,*)' medium-range intermolecular forces:',timel      
      write(*,*)' short-range intermolecular forces: ',times
	write(*,*)' change of subensemble:             ',time2
	write(*,*)' choice of neuigbours:              ',timev   
      write(*,*)' FURIR:                             ',timef      
      write(*,*)' ETERM:                             ',timee      
      write(*,*)' atom-atom intramolecular forces:   ',timen      
      write(*,*)' covalent bonds:                    ',timeb      
      write(*,*)' covalent angles:                   ',timea      
      write(*,*)' torsions:                          ',timet      
      write(*,*)' moving:                            ',timeg      
      write(*,*)'-------------------------------------------------'
      write(*,*)' Total time:                        ',time
	write(*,*)' node ',TASKID, NNN
	end if 
	write(*,*)
	write(*,*)'Acceptance ratio ',float(ICHU)/(float(ICHE)+0.001) 
	if(IWLI.lt.NWLI.or.NHIST.le.IHIST)write(*,*)
     + '!!! Expanded ensemble simulation not converged yet!!!!'
        write(*,*)'==================================================='
	write(*,*)' DISTRIBUTION OVER SUBENSEMBLES'
	write(*,'(2(a,I4),3x,f8.0,a)')' Series from ',IHIST,
     + ' to ',NHIST, SUMT,' config'
	write(*,*)' #  alpha     ita       prob.      bF     F(kJ/M)',
     +'   P-       P+'
	if(NDEL.eq.0)NDEL=1 
        if(NST.le.0)NST=1
        if(PLOW.le.1.d0/NST)PLOW=1.d0/NST
	do I=1,NE
	  PEE(I)=0.
	  do J=IHIST,NHIST
	    PEE(I)=PEE(I)+dfloat(ITE(I,J))/SUMT
	  end do
	end do
	do I=1,NE
          if(PEE(I).le.PLOW)then
            if(PEE(NE).le.PLOW)then
	      FEN(I)=(EE(I)-EE(NE))/NDEL                        ! Beta*F
            else
	      FEN(I)=(EE(I)-EE(NE)-dlog(PLOW/PEE(NE)))/NDEL     ! Beta*F
            end if
          else
            if(PEE(NE).le.PLOW)then
	      FEN(I)=(EE(I)-EE(NE)-dlog(PEE(I)/PLOW))/NDEL       ! Beta*F
            else
	      FEN(I)=(EE(I)-EE(NE)-dlog(PEE(I)/PEE(NE)))/NDEL    ! Beta*F
            end if
          end if
          FF=FEN(I)*0.001*ENERF/BETA                         ! F (kJ/M)
          IF(PEE(I).NE.0.) THEN
            if(I.gt.1) PTR0=IWALK(I,I-1)*NCEE/(PEE(I)*SUMT*IAVER)
            if(I.lt.NE)PTR1=IWALK(I,I+1)*NCEE/(PEE(I)*SUMT*IAVER)
          ELSE
            PTR0      = 0.
            PTR1      = 0.
          END IF
	  write(*,'(I3,f9.5,f9.4,f9.5,4f9.4)')
     +  I,EC(I),EE(I)/NDEL,PEE(I),FEN(I),FF,PTR0,PTR1
	end do
        write(*,*)'==================================================='
	write(*,*)' Distribution history:'
	do I=IHIST,NHIST
          if(HIST(3,I).eq.1.d0)write(*,'(a10,10x,80f8.3)')
     +  ' new n_ : ',(HIST(J,I),J=4,NE+3)
      write(*,'(i4,f8.2,f8.0,80f8.5)')
     +  I,HIST(1,I),HIST(2,I),(ITE(J,I)/HIST(2,I),J=1,NE)
	  TMP1(I)=ITE(1,I)/HIST(2,I)
	  TMP2(I)=ITE(NE,I)/HIST(2,I)
	end do                       
	call DISP(TMP1,PM1,PD1,IHIST,NHIST)
	call DISP(TMP2,PM2,PD2,IHIST,NHIST)
      write(*,*)'-----------------------------------------------'
	write(*,'(a,f10.7,a3,f10.7)')'  P(T) = ',PM1,'+/-',PD1
	write(*,'(a,f10.7,a3,f10.7)')'  P(0) = ',PM2,'+/-',PD2 
	if(PM1.gt.0.d0.and.PM2.gt.0.d0)then
	  FRM=(EE(NE)-EE(1)+dlog(PM1/PM2))/NDEL    ! -Beta*F
          FRM=-FRM*0.001*ENERF/BETA                ! F (kJ/M)
	  FRD=sqrt((PD1/PM1)**2+(PD2/PM2)**2)*0.001*ENERF/(BETA*NDEL)
  	  write(*,*)'==============================='
	  write(*,'(a,f10.4,a3,f10.4)')'  F    = ',FRM,'+/-',FRD 
  	  write(*,*)'==============================='
	end if
      IAV            = 2
      EBONDS	= 0.
      IC0=0
      STR(1:20)='#     alpha              '
      write(str(21:500),'(100x,100x,100x,100x,80x)')   
      DO ITYP        = 1,NTYPES
        NRBS           = NRB(ITYP)
        IF(NRBS.NE.0)  THEN
          ISHF           = ISADR(ITYP)
          NBBEG          = IADB(ITYP)
          NBEND          = NBBEG+NRBS-1
          IC0=17
          IAV0	= IAV+1
            DO I           = NBBEG,NBEND
            II             = IB(I)+ISHF
            JJ             = JB(I)+ISHF
            IAV            = IAV+2
            if(IC0.le.490)STR(IC0:IC0+5)=NM(II)(1:2)//'-'//NM(JJ)(1:2)
            IC0=IC0+17
          END DO! OF I
          if(IPRINT.ge.6)then
	      IF(IC0.GT.500)IC0=500
            write(*,*)
            write(*,*)' Distribution of bonds.     Type ',ITYP
	      write(*,*)STR(1:IC0)
	      do J=1,NE
	        write(*,3333)J,EC(J),(AWE(I,J),I=IAV0,IAV)
	      end do
	    end if                
      END IF
      END DO! OF ITYP
 3333 format(I5,f7.4,100(1x,f8.4,f8.3))
*
      write(str(21:500),'(100x,100x,100x,100x,80x)')   
      DO ITYP        = 1,NTYPES
      NRAS           = NRA(ITYP)
      IF(NRAS.NE.0)  THEN
      NABEG          = IADA(ITYP)
      NAEND          = NABEG+NRAS-1
      ISHF           = ISADR(ITYP)
      IC0=17
      IAV0	= IAV+1
      write(*,*)'-----------------------------------------------'
      DO I           = NABEG,NAEND
      II             = IA(I)+ISHF
      JJ             = JA(I)+ISHF
      KK             = KA(I)+ISHF
      IAV            = IAV+2
      if(IC0.le.485)STR(IC0:IC0+14)=NM(II)//'-'//NM(JJ)//'-'//NM(KK)
      IC0=IC0+17
      END DO! OF I
          if(IPRINT.ge.6)then
	      IF(IC0.GT.500)IC0=500
            write(*,*)
            write(*,*)' Distribution of angles'
	      write(*,*)STR(1:IC0)
	      do J=1,NE
	        write(*,3334)J,EC(J),(AWE(I,J),I=IAV0,IAV)
	      end do
	    end if                
 	END IF
      END DO! OF ITYP
 3334 format(I5,f7.4,100(1x,f8.2,f8.3))
*
	write(str(21:500),'(100x,100x,100x,100x,80x)')  
	DO ITYP        = 1,NTYPES
      NRTS           = NRT(ITYP)
      IF(NRTS.NE.0)  THEN
      NTBEG          = IADT(ITYP)
      NTEND          = NTBEG+NRTS-1
      ISHF           = ISADR(ITYP)
      IC0=17
      IAV0	= IAV+1
          write(*,*)'-----------------------------------------------'
      DO I           = NTBEG,NTEND
      II             = IT(I)+ISHF
      JJ             = JT(I)+ISHF
      KK             = KT(I)+ISHF
      LL             = LT(I)+ISHF
      IAV            = IAV+2
      if(IC0.le.485)STR(IC0:IC0+14)=NM(II)//'-'//NM(JJ)//'-'//NM(LL)
      IC0=IC0+17
      END DO! OF I
          if(IPRINT.ge.6)then
	      IF(IC0.GT.500)IC0=500
            write(*,*)
            write(*,*)' Distribution of torsions'
	      write(*,*)STR(1:IC0)
	      do J=1,NE
	        write(*,3333)J,EC(J),(AWE(I,J),I=IAV0,IAV)
	      end do
	    end if                
 	END IF
      END DO! OF ITYP
*
      write(str(21:500),'(100x,100x,100x,100x,80x)')   
      IC0=17
      IAV0	= IAV+1
      if(IPRINT.ge.7)write(*,*)'------------------------------------'
      DO I           = 1,NTYPES
        IAV            = IAV+2
        if(IPRINT.ge.7)then
          PRINT "('*** SOME ENERGY CONTRIBUTIONS: *** - TYPE:',I2)",I
	  write(STR(IC0:IC0+1),'(I2)')I
          if(IC0.le.485)STR(IC0+2:IC0+17)=' : PES14   PLJ14'
          IC0=IC0+17
        end if
      END DO
      IF(IC0.GT.500)IC0=500
      if(IPRINT.ge.7)then
        write(*,*)
        write(*,*)' Distribution of some intramolecular energies'
	  write(*,*)STR(1:IC0)
	  do J=1,NE
	    write(*,3335)J,EC(J),(AWE(I,J),I=IAV0,IAV)
	  end do
	end if
 3335 format(I5,f7.4,50(1x,3f8.3))
*
      write(str(21:500),'(100x,100x,100x,100x,80x)')  
      IC0=17
      IAV0	= IAV+1 
      if(IPRINT.ge.5)then
        write(*,*)
        write(*,*)'-----------------------------------------------'
      end if
      DO ITYP        =    1,NTYPES
        DO JTYP        = ITYP,NTYPES
          I              =  MDX(ITYP,JTYP)
          IAV            = IAV+2
          if(IC0.le.485)STR(IC0:IC0+12)=NAME(ITYP)//'-'//NAME(JTYP)
          IC0=IC0+17
        END DO! OF JTYP
      END DO! OF ITYP
      IF(IC0.GT.500)IC0=500
      if(IPRINT.ge.5)then
        write(*,*)
        write(*,*)' Distribution of intermolecular energies'
	  write(*,*)STR(1:IC0)
	  NP=NTYPES*(NTYPES+1)/2
	  do IP=1,NP
	    TINN(IP)=0.
	  end do
	  do J=1,NE
	    write(*,3336)J,EC(J),(AWE(I,J),I=IAV0,IAV)
	  end do           
	end if                
 3336 format(I5,f7.4,50(1x,2f8.2))
*
      write(str(19:500),'(100x,100x,100x,100x,82x)')   
	IAV0	= IAV+1
      write(*,*)'-----------------------------------------------'
	write(*,*)
      IAV            = IAV+6
	STR(13:19)=' TEMP  '
	STR(20:28)=' EPinter '
	STR(29:37)=' EPintra '
	STR(38:46)=' EPtotal '
	STR(47:55)=' PRES    '
	STR(56:64)=' Esolv   '
	STR(65:73)=' TS      '
      write(*,*)
      write(*,*)' Distribution of global thermodynamic properties'
	write(*,*)STR(1:80)
	ET0=AWE(IAV0+1,NE)+AWE(IAV0+2,NE)
        if(LSHEJK)then
           IAVP=IAV0+4
        else
           IAVP=IAV0+3
        end if
	do J=1,NE                                                   
	  ETOT=AWE(IAV0+1,J)+AWE(IAV0+2,J)
	  ESOLV=NOP*(ETOT-ET0)/NDEL
          FF=FEN(J)*0.001*ENERF/BETA                         ! F (kJ/M)
          SS=ESOLV-FF
	  write(*,'(I3,f8.4,f8.2,3f9.3,5f9.2)')
     +  J,EC(J),AWE(IAV0,J),AWE(IAV0+1,J),AWE(IAV0+2,J),
     +  ETOT,AWE(IAVP,J),ESOLV,SS  
	end do
	PP1=AWE(IAVP,1)
	PPM=AWE(IAVP,NE)
	if(LNPT)then
          VOL1            = AWE(IAV,1)
          VOLM            = AWE(IAV,NE)
	  PVTERM1	= -TRPRES*VOL1*AVSNO*1.d-28/NOP
	  PVTERM2	= -TRPRES*VOLM*AVSNO*1.d-28/(NOP-NDEL)
	else
	  PVTERM1	= -PP1*VOL*AVSNO*1.d-28/NOP
	  PVTERM2	= -PPM*VOL*AVSNO*1.d-28/(NOP-NDEL)
	end if 
	PVTERM=0.5*(PVTERM1+PVTERM2)
*	write(*,*)' Pressure= ',PRESA, 'atm'
*	if(.NOT.LNPT)then
          PKT=0.001*ENERF/BETA
          PKT2=0.001*ENERF*NOP/(BETA*(NOP-NDEL))
*	  write(*,*)' PV-term 1 = ',PVTERM1,' kJ/M' 
*	  write(*,*)' PV-term 2 = ',PVTERM2,' kJ/M' 
	  write(*,*)' kT-term   = ',PKT,' kJ/M'
*	  write(*,*)' kT-term 2 = ',PKT2,' kJ/M'
	  write(*,671)' Free energy:',FRM,'+/-',FRD,' kJ/M' 
*	  write(*,671)' Alt. Fr.en.:',FRM+PVTERM2+PKT2,'+/-',FRD,' kJ/M' 
*	else
*	  PKT=-0.001*ENERF*NDEL*dlog(TPRES*VOL*BETA/(UNITP*NDEL))/BETA
*	  write(*,*)' kT-term = ',PKT,' kJ/M'
*	  write(*,671)' Free energy:',FRM+PKT,'+/-',FRD,' kJ/M'
*	end if 
 671	format(a,f12.4,a,f9.4,a)
	write(*,*)'----------------------------------------------------'
      write(*,'(17x,3f9.3,4f9.2)')AV1,AV2,AV3,AV4,AV5,AV6
	if(IPRINT.lt.4)return
*
      write(*,*)
      write(*,*)' Fractional temperatures'
      STR(14:80)=
     +'   Dens     boxX    boxY   boxZ    Ttrans    Trot     Tint'    
	write(*,*)STR(1:75)
	do J=1,NE
          AVX = AWE(IAV,J)
          AVY = AWE(IAV+1,J)
          AVZ = AWE(IAV+2,J)
          AV7            = AVX*AVY*AVZ
*  Mass of deleted particles
	  TOTMEE = 0
	  do ITYP=1,NTYPES
	    if(IEE(ITYP).eq.0)
     +  TOTMEE = TOTMEE + SUMMAS(ITYP)*NSPEC(ITYP)/(AVSNO*1000.)
	  end do
          DENS	= (TOTMAS-(1.-EC(J))*TOTMEE)*1.d-3/(AV7*UNITL**3)
	  write(*,'(I5,f8.4,f9.4,3f9.2,3f9.3)')
     +  J,EC(J),DENS,
     +  AVX,AVY,AVZ,(AWE(I,J),I=IAV+3,IAV+5)
	end do
      write(*,*)'----------------------------------------------------'
*     
      write(*,*)' Table of transitions'
      write(*,'(4x,50i4)')(I,I=1,NE)
      do J=1,NE
        write(*,'(51i4)')J,(IWALK(J,I),I=1,NE)
      end do          
      RETURN
      END 

