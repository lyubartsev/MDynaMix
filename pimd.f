*     MDynaMix v. > 5.1
*
*     File pimd.f
*     ------------
*
C     This file contains subroutines relevant to Path Integral Molecular Dynamics
*
*====================================================================
*
*   1. PIMD - Integration of motion in Path Integral Molecular Dynamics
*
*   Centroid MD in the adiabatic approximation and
*   use of double time step algorithm
*   (will be described in more details in some future)
*
*   This subroutine repeats in many features the DOUBLE subroutine in mdstep.f file
C  
      SUBROUTINE PIMD
      include "common_inc.f"
*
*   1.1 Local definitions
*
C   LRDF defines wheter to calculate RDFs
      LOGICAL LRDF,DONE  
C   this is temporary array to keep long-range contributions to virial
      real*8 TRLL(3)
C   RDF should not be calculated during initialisation stage (2.2)
*      write(*,*)' node ',JBEAD,' enter PIMD'
      LRDF=.false.
      NFREQ2=NFREQ/2+1
      VIRA(1)=0.
      VIRA(2)=0.
      VIRA(3)=0.
      NEKS=0
      QEKS=0.
      TEMPNOD=0.
      IPIMA=0
      do I=1,NSITES
         PIRS(I)=0.
      end do
      do I=1,NSTOT
        do J=1,JNH
          THERMS(I,J)=0.
        end do
      end do
*
*   1.2 Preparation
*   ---------------
*   
*   1.2.1  Compose list of neigbours       
      call CHCNB(LRDF)      ! mdstep.f 
*   1.2.2  Calculate centre-of-masses
      CALL GETCOM    
*   1.2.3  Calculate all forces (should be done before entering the cycle)
      LRDF=.true.  
C   slow forces  (important: slow forces calculated before fast ones!)
      call SLFORCES
C   Fast forces
      call FFORCES(.false.)
C   Scale down the  forces
      do I=1,NSTOT
        GX(I)=GX(I)/NBEADS
        GY(I)=GY(I)/NBEADS
        GZ(I)=GZ(I)/NBEADS
        HX(I)=HX(I)/NBEADS
        HY(I)=HY(I)/NBEADS
        HZ(I)=HZ(I)/NBEADS
        TEMBS(I)=0.
      end do
C  Forces between PI beads
      call PIFORCES
*      write(*,'(i3,a8,7e12.5)')JBEAD,' PIRR = ',(PIRR(IS),IS=1,NSITES)
*   1.2.4  Calculate different contributions to pressure and energy
C     (this is done mostly for control purposes)
      timeg0	= cputime(0.d0)
      if(IPRINT.ge.7)then
	PELL=PELS1+PELS2+SPE
        PENL=(PE1+PE2+PELL)*PERMOL
        PEL=PELL*PERMOL
        PENS=PINT*PERMOL
        PEKI=TKE*PERMOL
        ETOT=PENL+PENS+PEKI
        WRITE(6,'(4(a6,f12.5))')
     +' Eex= ',PENL,' Ein= ',PENS,' Eel= ',PEL,' Etot=',ETOT 
      end if
*
*   1.3  Double time step algorithm
*   -------------------------------
*   1.3.1 Organize the cycle (till 2.4.9)
      NSTEP=0
 1    NSTEP=NSTEP+1
      MSTEP=NSTTOT+NSTEP
*   1.3.2 First thermostat integration: to half of the long time step
      call PI_SCALING(1)
*   1.3.3 First integration of slow forces
C     Velocities are corrected as due to slow forces acting during 
C     half of the long time step
      DO I           =  NAB(TASKID),NAE(TASKID)
	ITYP	= ITYPE(I)
C   LMOVE define which molecules can move (normally, LMOVE=.true.)
	if(LMOVE(ITYP))then 
          VX(I)          =   VX(I)+HTSTEP*GX(I)
          VY(I)          =   VY(I)+HTSTEP*GY(I)
          VZ(I)          =   VZ(I)+HTSTEP*GZ(I)
	end if
      END DO! OF I
*
*    1.3.4  Integration of the fast forces
*    ------------------------------------- 
*    1.3.4.1  Cycle over short time steps
C  This set zeros for current averages bond/angles and their energies
C  These averages are taken over the molecules and over short time
C  steps within this large step 
      call ZEROAV
      QEKS=0
      DO NSTEB       =  1,NFREQ
        DO I           =  NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then 
*    1.3.4.2  Velocities at t(short)+1/2
            VX(I) =   VX(I)+HTSTEB*HX(I)         ! V(t+1/2)
            VY(I) =   VY(I)+HTSTEB*HY(I)
            VZ(I) =   VZ(I)+HTSTEB*HZ(I)
          end if
        end do
*   1.3.4.3. Computation of centroid momentum  FX,FY,FZ
*        call CMD_MOM                     ! mpi communication
*   1.3.4.4  Nose themostat for local motions of beads
        call BEADNOSE
*    1.3.4.5  Coordinates at t(short)+1
        DO I           =  NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then 
            SX(I) =  VX(I)*TSTEB*MASSDI(I)+SX(I)     ! X(t+1)
            SY(I) =  VY(I)*TSTEB*MASSDI(I)+SY(I)
            SZ(I) =  VZ(I)*TSTEB*MASSDI(I)+SZ(I)
	  end if
        END DO! OF I
        timeg	= timeg+cputime(timeg0)   
*  1.3.4.6  Recalculate centre-of-mass coordinates
	CALL GETCOM         ! services.f
*
*  1.3.5 Recalculate forces
*  ------------------------
*  1.3.5.1 After "ICHNB" long time steps - recalculate list of neighbours
C  This is done in the middle of long time step  
	if(NSTEB.eq.NFREQ2.and.mod(MSTEP,ICHNB).eq.0)call CHCNB(LRDF)     
	if(NSTEB.eq.NFREQ)then
*  1.3.5.2 Each long time step - calculate slow forces
          CALL SLFORCES
*  1.3.5.3 Claculate fast forces
C         The logical parameter signals whether to calculate averages 
          call FFORCES(.true.)
          do I=1,NSTOT
            GX(I)=GX(I)/NBEADS
            GY(I)=GY(I)/NBEADS
            GZ(I)=GZ(I)/NBEADS
            HX(I)=HX(I)/NBEADS
            HY(I)=HY(I)/NBEADS
            HZ(I)=HZ(I)/NBEADS
          end do
C  Forces between PI beads
          call PIFORCES
        else
          call FFORCES(.false.)
          do I=1,NSTOT
            HX(I)=HX(I)/NBEADS
            HY(I)=HY(I)/NBEADS
            HZ(I)=HZ(I)/NBEADS
          end do
          call PIFORCES
        end if
        QEKS = QEKS + QEK
*   
*   1.3.6 Conclude integration of fast forces
*   -----------------------------------------
C   (continued in the beginning of the cycle, see p.2.3.2    
        timeg0	= cputime(0.d0)
	if(IPRINT.ge.8)write(*,'(a,I3,f12.3)')' 1.st.',NSTEB,TEMP
        DO I           =  NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then 
            VX(I)          =  VX(I)+HTSTEB*HX(I)     !  V(t+1)
            VY(I)          =  VY(I)+HTSTEB*HY(I)
            VZ(I)          =  VZ(I)+HTSTEB*HZ(I)
          end if
	END DO! OF I      
        call GETEMP
	if(IPRINT.ge.8)write(*,'(a,2I3,f12.3)')' sh.st.',JBEAD,NSTEB,TEMP
*
      END DO! OF NSTEB 
*
*   1.3.7 Conclude integration of slow forces
*   -----------------------------------------  
*   1.3.7.1 Velocities at t(long)+1
        DO I           = NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then
            VX(I)          =  VX(I)+HTSTEP*GX(I)
            VY(I)          =  VY(I)+HTSTEP*GY(I)
            VZ(I)          =  VZ(I)+HTSTEP*GZ(I)
	  end if
        END DO! OF I
      timeg	= timeg+cputime(timeg0)
*   1.3.7.2 Conclude integration of NHC
      call PI_SCALING(2)
*   1.3.6.3 Report velocities to the master node
      timeg0	= cputime(0.d0)
*
*   1.4 calculate averages for intermediate output
*   ----------------------------------------------
* remove excess COM momenta if specified
      if(ICMOM.gt.0.and.mod(MSTEP,ICMOM).eq.0)call CHKMOM
*   1.4.1 Different contribution to the kinetic energy 
        CALL GETKIN       ! services.f
*  Accumullation of data from the nodes
*   1.4.4 Acccumulate averages
        TEMPNOD=TEMPNOD+TEMP
        call PI_AVER
        PEST = PELS1+PELS2+SPE
        PE = PE1+PE2+PEST
        if(mod(MSTEP,IAVER).eq.0)CALL GETAVR(1)  ! aver.f
*   1.4.5 Dump trajectory
	if(mod(MSTEP,NTREK).eq.0)CALL TRACE      ! restart.f
*   1.4.6 Calculate different contributions to pressure and energy
        PENL = PE*PERMOL                ! intermolecular energy
        PEL= PEST*PERMOL                ! electrostatic energy          
        PENS=PINT*PERMOL                ! intramolecular energy
        PEKI=TKE*PERMOL                 ! classical kinetic energy
        POTEN=PENL+PENS                 ! total potential energy
        ETOT =POTEN+PEKI                ! total energy   
        WIRSUM=WIRS+WIRSS
        TRYCKE	= PEL*UNITP/(3.*VOL*PERMOL)     ! electrostatic
        TRYCKL	= (VIR1+VIR2)*UNITP/(3.*VOL)    ! LJ
        TRYCKB	= VIRB*UNITP/(3.*VOL)           ! bonds
        TRYCKS	= WIRSUM*UNITP/(3.*VOL)         ! molecular
        TRYCKP	= VIRPI*UNITP/(3.*VOL)         ! PI springs
*        TRYCKA	= (VIRA(1)+VIRA(2)+VIRA(3))*UNITP/(3.*VOL) ! molecular
        TRYCKD	= VIRD*UNITP/(3.*VOL)           ! long-range corr.
        DENS	= TOTMAS*1.d-3/(VOL*UNITL**3)
C  Quantum kinetic energy
        QEK = QEKS / NFREQ
        TKEQ = 1.5*NSTOT*NBEADS*BOLTZ*TRTEMP/UNITE-QEK
*   1.4.7 Current output
C OUTPUT for visualisation
        if(LVISUAL)then
          if(mod(MSTEP,IOUTPUT).eq.0)then
            write(*,'(a)')'@mm coord.start'
            do I=1,NSTOT
              write(*,'(a,3f12.4)')'@mm ',SX(I),SY(I),SZ(I)
            end do
            write(*,'(a)')'@mm coord.end'
          end if
          WRITE(*,'(a,I7,4(1X,F9.3),2(1X,f9.2),f9.4,3f10.4)') 
     +    '@mm E ',MSTEP,POTEN,PENS,PEL,ETOT,TEMP,TRYCKM,DENS,
     +     BOXL,BOYL,BOZL
        end if   ! if(LVISUAL
        if(mod(MSTEP,IOUTPUT).eq.0)then
          if(IPRINT.ge.5)then
            WRITE(6,'(I8,4(1X,F9.3),2(1X,f9.2),f9.4)')
     +      MSTEP,POTEN,PENS,PEL,ETOT,TEMP,TRYCK,DENS
          if(JBEAD.eq.0)write(*,'(a,f10.3,a,f10.3)')
     +    'Eharm: ',QEK*PERMOL,'  Ekin_Q:',TKEQ*PERMOL
        end if
        if(IPRINT.ge.6)then
          write(*,'(2f7.4,f9.2,4f7.0,I7,I9)')
     + SC,SCL,TRYCKM,TTR,TROT,TINT,TEMPV,MBSH,MBLN   
          if(JBEAD.eq.0)write(*,'(a,80f8.4)')
     +    'RRav: ',(PIRR(I),I=1,NSITES)
   	  if(LSEP)write(*,'(4f11.1,3f11.3)')
     + TRYCKX,TRYCKY,TRYCKZ,(TRYCKX+TRYCKY+TRYCKZ)/3.,BOXL,BOYL,BOZL
          if(IEXT.eq.1)write(*,'(2(a,e14.6))')
     +    ' E abs: ',(EABS+EABSS)*PERMOL,' Ext.f. ',EFIELD
        end if
 	  if(IPRINT.ge.7.and.mod(MSTEP,IOUTPUT).eq.0)
     +  WRITE(*,'(8f10.1)')TRYCKL,TRYCKE,TRYCKB,TRYCKS,TRYCKP,TRID
          if(IPRINT.ge.6.and.LSCTD)
     +    write(*,'(8f10.5)')(TEMPR(I),SCM(I),I=1,NTYPES)
*** Put you own subroutine for evaluation of whatever you want here:
***     call USER
***
        end if  ! if(mod(MSTEP
*   1.4.8  Intermediate output after "NAVER" steps   
      if(mod(MSTEP,NAVER).eq.0.and.JBEAD.eq.0)CALL GETAVR(2)
*   1.4.9  Dump restart file   
C   Important - all nodes
      if(LDMP.and.mod(MSTEP,NDUMP).eq.0)call RESTRT(2)
      timeg	= timeg+cputime(timeg0)
*
*   1.4.10 Conclude MD step
*   Check STOP file (visual regime only)
      if(LVISUAL)then
        open(unit=77,file='MD_STOP',status='old',err=10)
        close(77)
        return
      end if
 10   if(NSTEP.lt.NSTEPS)go to 1
*<---------------------
*   PIMD termostate analysis
      write(*,*)' Node ',JBEAD,' averager temperature ',TEMPNOD/NSTEPS
      if(IPRINT.ge.5)then
        write(*,*)' particle temperatures and thermostats factors:'
        do I=1,NSTOT
          write(*,'(i6,f11.4,10f10.6)')
     +    I,TEMBS(I)/IPIMA,(THERMS(I,J)/(FNH*IPIMA),J=1,JNH)
        end do
        write(*,*)' average lengths of springs:'
        do J=1,NSITES
          write(*,'(a,i4,f12.4)')
     +  ' site ',J,PIRS(J)/IPIMA
        end do
      end if
      RETURN
      END
C
C   2. Calculate PI harminic forces between the beads
C
*
*==================== PIFORCES =============================================
*
      subroutine PIFORCES
      include "common_inc.f"
      include "mpif.h"
      real*8 PIRRT(NS)
      VIRPI=0.
C     
* 2.1  Retrive coordinates from upper bead (copied into array OX,,,) 
      call PI_SENDCOORD(1)
      do I=1,NSITES
         PIRRT(I)=0.
      end do
      QEK=0.
C     
* 2.2 step forward
      do I=1,NSTOT
        IS = NSITE(I) 
        BX        = SX(I)-OX(I)
        BY        = SY(I)-OY(I)
        BZ        = SZ(I)-OZ(I)
        call PBC(BX,BY,BZ)
        BSQ       = BX*BX+BY*BY+BZ*BZ
        BBB        = DSQRT(BSQ)
        R1I       = 1.D0/BBB
C    Potential parameter
        DBND      = BBB*BEADFAC*MASS(IS)
C    Energy of the bond
        EBDD	    = DBND*BBB
C    Force (abs. value)
        FORB      =  -2.0D0*DBND*R1I
        HX(I)     = HX(I)+BX*FORB
        HY(I)     = HY(I)+BY*FORB
        HZ(I)     = HZ(I)+BZ*FORB
        QEK = QEK + EBDD            ! PI kin energy
        PIRRT(IS)=PIRRT(IS)+BBB       ! Average spring length
        VIRPI = VIRPI + FORB*BSQ
	VIRFX = VIRFX+FORB*BX**2*NBEADS    ! will be scaled back by 1/NBEADS  
	VIRFY = VIRFY+FORB*BY**2*NBEADS  
	VIRFZ = VIRFZ+FORB*BZ**2*NBEADS
      end do
*  2.3. Retrive coordinates from upper bead (copied into array OX,,,) 
      call PI_SENDCOORD(-1)
C     
      do I=1,NSTOT
        IS = NSITE(I) 
        BX        = SX(I)-OX(I)
        BY        = SY(I)-OY(I)
        BZ        = SZ(I)-OZ(I)
        call PBC(BX,BY,BZ)
        BSQ       = BX*BX+BY*BY+BZ*BZ
        BBB        = DSQRT(BSQ)
        R1I       = 1.D0/BBB
C    Potential parameter
        DBND      = BBB*BEADFAC*MASS(IS)
C    Energy of the bond
        EBDD	    = DBND*BBB
C    Force (abs. value)
        FORB      =  -2.0D0*DBND*R1I
        HX(I)     = HX(I)+BX*FORB
        HY(I)     = HY(I)+BY*FORB
        HZ(I)     = HZ(I)+BZ*FORB
        QEK = QEK + EBDD          ! PI kin energy
      end do
      NSIT1=NSITES+1
      PIRRT(NSIT1)=QEK
      call MPI_ALLREDUCE(PIRRT,PIRR,NSIT1,MPI_DOUBLE_PRECISION,
     +     MPI_SUM,MPI_COMM_WORLD,ierr)
      QEK=0.5*PIRR(NSIT1)
      do IS=1,NSITES
        ITYP=ITS(IS)
        PIRR(IS)=PIRR(IS)/(NSPEC(ITYP)*NBEADS)
      end do
*      write(*,*)' Eharm',QEK*PERMOL,JBEAD
      return
      end
C
C   3. Nose thermostate for beads
C   Note: separate thermostat for each atom !!!
C       todo: update by the chain of thermostat
*
*==================== BEADNOSE ===============================
*
      subroutine BEADNOSE
      include "common_inc.f"
      include "mpif.h"
      real*8 FEK(NTOT),BUF(NTOT)
*  3.1 Compute faked kinetic energy
      do I=1,NSTOT
        BUF(I)=0.5*(VX(I)**2+VY(I)**2+VZ(I)**2)*MASSDI(I)
      end do
      call MPI_ALLREDUCE(BUF,FEK,NSTOT,MPI_DOUBLE_PRECISION,
     +     MPI_SUM,MPI_COMM_WORLD,ierr)
      do I=1,NSTOT 
        ATEMP=FEK(I)*CONVEQ
        DKE = ATEMP/TRTEMP-1.d0
        THERMO(I,1)=THERMO(I,1) + DKE*DQQT*TSTEB
        THERMS(I,1)=THERMS(I,1)+THERMO(I,1)**2
C  Nose-Hoover chain
        if(JNH.gt.1)then
          do J=2,JNH
            THERMO(I,J) = THERMO(I,J) + (THERMO(I,J-1)**2-FNH)*TSTEB
            THERMO(I,J-1) = THERMO(I,J-1)*(1.-THERMO(I,J)*TSTEB)
            THERMS(I,J)=THERMS(I,J)+THERMO(I,J)**2
          end do          
        end if
        SCV=1.-THERMO(I,1)*TSTEB 
        if(IPRINT.ge.8)
     +    write(*,*)' PITH:',I,ATEMP,(THERMO(I,J),J=1,JNH),SCV
C        VX(I)=(VX(I)-FX(I))*SCV+FX(I)
C        VY(I)=(VY(I)-FY(I))*SCV+FY(I)
C        VZ(I)=(VZ(I)-FZ(I))*SCV+FZ(I)
        VX(I)=VX(I)*SCV
        VY(I)=VY(I)*SCV
        VZ(I)=VZ(I)*SCV
        TEMBS(I)=TEMBS(I)+ATEMP
      end do
      do I=1,NSITES
        PIRS(I)=PIRS(I)+PIRR(I)
      end do
      IPIMA=IPIMA+1
      return
      end
C
C   4. Get coordinates of neighbouring node
C   ------------------------------------
*
*==================== PI_SENDCOORD ===============================
*
      subroutine PI_SENDCOORD(IDIR)
      include "common_inc.f"
      include "mpif.h"
      real*4 BUFF(NBUFF),BUF2(NBUFF)
      integer ISTS(MPI_STATUS_SIZE)
      if(IDIR.gt.0)then
        ISEND=JBEAD+1
        if(ISEND.ge.NBEADS)ISEND=0
        IRECV=JBEAD-1
        if(IRECV.lt.0)IRECV=NBEADS-1
      else
        ISEND=JBEAD-1
        if(ISEND.lt.0)ISEND=NBEADS-1
        IRECV=JBEAD+1
        if(IRECV.ge.NBEADS)IRECV=0
      end if
      NST3=3*NSTOT
      IMA=0
      do I=1,NSTOT
        BUFF(I)=SX(I)
        BUFF(I+NSTOT)=SY(I)
        BUFF(I+2*NSTOT)=SZ(I)
      end do
*      write(*,*)' node ',JBEAD,' in SENDCRD'
      if(mod(JBEAD,2).eq.0)then
        call MPI_SEND(BUFF,NST3,MPI_REAL,ISEND,IMA,MPI_COMM_WORLD,ierr)
        call MPI_RECV(BUF2,NST3,MPI_REAL,IRECV,IMA,
     +  MPI_COMM_WORLD,ISTS,ierr)
      else
        call MPI_RECV(BUF2,NST3,MPI_REAL,IRECV,IMA,
     +  MPI_COMM_WORLD,ISTS,ierr)
        call MPI_SEND(BUFF,NST3,MPI_REAL,ISEND,IMA,MPI_COMM_WORLD,ierr)
      end if
*      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      do I=1,NSTOT
        OX(I)=BUF2(I)
        OY(I)=BUF2(I+NSTOT)
        OZ(I)=BUF2(I+2*NSTOT)
      end do
      return
      end      
C
C   This compute molecular center of mass of centroids
C
*
*=================== PI_COM ======================================
*
      subroutine PI_COM(XPC,YPC,ZPC)
      include "common_inc.f"
      include "mpif.h"
      dimension XPC(NPART),YPC(NPART),ZPC(NPART)
      double precision BUFF(NBUFF),BUF2(NBUFF)
      do I=1,NOP
        BUFF(I)=X(I)
        BUFF(I+NOP)=Y(I)
        BUFF(I+2*NOP)=Z(I)
      end do
      IBUF=3*NOP
      call MPI_ALLREDUCE(BUFF,BUF2,IBUF,MPI_DOUBLE_PRECISION,
     +     MPI_SUM,MPI_COMM_WORLD,ierr)
      do I=1,NOP
         XPC(I)=BUF2(I)/NBEADS
         YPC(I)=BUF2(I+NOP)/NBEADS
         ZPC(I)=BUF2(I+2*NOP)/NBEADS
      end do
      return
      end
C
C   This compute centroid momentum by summation momentum from all nodes
C
*
*=================== CMD_MOM ======================================
*
      subroutine CMD_MOM
      include "common_inc.f"
      include "mpif.h"
      real*4 BUFF(3*NTOT),BUF2(3*NTOT)
      do I=1,NSTOT
        BUFF(I)=VX(I)
        BUFF(I+NSTOT)=VY(I)
        BUFF(I+2*NSTOT)=VZ(I)
      end do
      IBUF=3*NSTOT
      call MPI_ALLREDUCE(BUFF,BUF2,IBUF,MPI_REAL,
     +     MPI_SUM,MPI_COMM_WORLD,ierr)
      do I=1,NSTOT
         FX(I)=BUF2(I)
         FY(I)=BUF2(I+NSTOT)
         FZ(I)=BUF2(I+2*NSTOT)
      end do
      return
      end
*
*=============== SCALING ==============================================
*
      SUBROUTINE PI_SCALING(IAL)
*
*  Temperature / pressure scaling
*  Path Integral version
*  ----------------------------------------
C
C  IAL = 1 - first half-step in double time step algorithm
C  IAL = 2 - second half-step in double time step algorithm
C
*
*  1 Definitions
*  ---------------
      include "common_inc.f"   
      parameter (FP3=0.5/3.)
      real*8 DPE(6),DPES(6),DTE(NTPS)
      data DKEO/0./,ISKK/0/
*
*  3 Recalculate temperatures and pressures
*  ---------------
*  3.3.1 Calculate temperature    
      call GETEMP
*  3.3.2 Calculate pressure  
C        This collect virial for "atomic" pressure
      DPE(1) = (VIR1+VIR2+VIRB+PELS1+PELS2+SPE+VIRD+
     +            VIRA(1)+VIRA(2)+VIRA(3))/NBEADS
C        This collect virial for "molecular" pressure
      DPE(2) = (VIR1+VIR2+PELS1+PELS2+SPE+VIRD+WIRS+WIRSS)/NBEADS
CD      write(*,'(5e13.5)')VIR1,VIR2,PELS1,PELS2,SPE,VIRD,WIRS,WIRSS
      DPE(3) = VIRPI 
      DPE(4) = VIRFX/NBEADS
      DPE(5) = VIRFY/NBEADS
      DPE(6) = VIRFZ/NBEADS
      call ALLSUM(DPE,DPES,6)
      VIRPI =  DPES(3)
      VIRSUM = DPES(1)+VIRPI
      WIRSUM = DPES(2)+VIRPI
      VIRFX = DPES(4)
      VIRFY = DPES(5)
      VIRFZ = DPES(6)
C        Kinetic (ideal gas) contributions to pressure
      TRID	= 2.*TKE*UNITP*NBEADS/(3.*VOL)
      TRIDM	= NOP*BOLTZ*TEMP*1.d-5*NBEADS/(VOL*UNITL**3)
C        Pressure in atm and its projections
      TRYCK	= VIRSUM*UNITP/(3.*VOL)+TRID
*      TRYCKM	= WIRSUM*UNITP/(3.*VOL)+TRIDM
      TRYCKM=0.
C        Projections are calculated for "atomic" pressure
      TRYCKX	= (VIRX+VIRFX+VIRD/3.+VIRA(1))*UNITP/VOL+TRID
      TRYCKY	= (VIRY+VIRFY+VIRD/3.+VIRA(2))*UNITP/VOL+TRID
      TRYCKZ	= (VIRZ+VIRFZ+VIRD/3.+VIRA(3))*UNITP/VOL+TRID
CD	if(TASKID.eq.MAST)write(*,'(5f13.2)')
CD     +TRYCK,TRYCKX,TRYCKY,TRYCKZ,(TRYCKX+TRYCKY+TRYCKZ)/3.
C
C     Now temperature TEMP and pressure TRYCK are defined for the 
C     current configuration, as well as their components
*
*   3.3.3 Calculate deviation from thermostat T and P
      DKE         = TEMP/TRTEMP-1.d0               !  total temperature
      do I=1,NTYPES
CD         write(*,*)' temp ',I,TEMPR(I),TASKID
         DTE(I)   = TEMPR(I)/TRTEMP-1.d0           !  for each species
      end do
C    in the case of constrained dynamics, pressure is defined by 
C    "molecule" algorithm - it corresponds scaling of molecular COM
C    for flexible molecules, scaling of atom positions is employed,
C    which corresponds to the "atomic" pressure 
C    Exception: case of separate pressure control in each direction 
C    (LSEP=.true.). Then pressure in each direction is determined from
C    "atomic" pressure. 
      if(LSHEJK)then
	 DPE(1) = TRYCKM-TRPRES
      else
         DPE(1) = TRYCK -TRPRES
      end if
      DPE(2) = TRYCKX-TRPRES
      DPE(3) = TRYCKY-TRPRES
      if(LHEX)then
        DPE(2)=0.5*(DPE(2)+DPE(3))
        DPE(3)=DPE(2)
      end if
      DPE(4) = TRYCKZ-TRPRES
*
*   3.4  First integration of Nose-Hoover equations  (t->t+1/2)
*   -----------------------------------------------------------
      if(IAL.eq.1)then
        if(LNVT)then
*   3.4.1 Correction ksi due to thermostat
C         coefficients DQT,DQP and RTP were defined in main.f, p.6.2
C
C    Nose eqn is:    dP/dt   -> F - P*ksi/Q
C                    dksi/dt -> 2*Ekin - Nf*kT
C    Here ksi/Q = SC
C             Q = Nf*kT*tau**2  ; tau -> QT (i.u.)
C                                 Q/(Nf*kT) -> 1./DQT 
C    So     d(SC)/DT = (2*Ekin/Nf*kT-1)*DQT = DKE*DQT
C
C    This is separate thermostat for each species 
          if(LSCTD)then
            SCA = 0.
            do I=1,NTYPES
              SCM(I) = SCM(I) + DTE(I)*DQT*HTSTEP     !  DQT -thermostat mass
              SCA = SCA + SCM(I)*EKIN(I) 
            end do
            SC = SCA/TKE
C   common thermostat
          else
 	    SC = SC + DKE*DQT*HTSTEP
          end if
          if(LNPT)then
*   3.4.2 Barostat corrections (NPT)
*   -------------------------------
*   3.4,2,1 Correction ksi due to barostat
C
C   Additional correction to ksi:
C       d(ksi)=ita**2/W-kT
C       ita/W -> SCL
C     so d(ksi)  ->   SCL**2*(W/Q) - kT/Q 
C                       
            if(LSCTD)then      
              SCA = 0.
              do I=1,NTYPES
                SCM(I) = SCM(I) + (SCL**2*RTP-DQT/FNST)*HTSTEP
                SCA = SCA + SCM(I)*EKIN(I) 
              end do
              SC = SCA/TKE
            else
              SC = SC + (SCL**2*RTP-DQT/FNST)*HTSTEP    !  RTP=DQT/DQP
            end if
*   3.4.2.2 Correction ita 
C         TKE is kinetic energy
            if(LSEP)then
              SCLX = SCLX+DQP*(VOL*DPE(2)+2.*TKE/FNST-SCLX*SC)*HTSTEP
              SCLY = SCLY+DQP*(VOL*DPE(3)+2.*TKE/FNST-SCLY*SC)*HTSTEP
              if(LHEX)then
                 SCLX=0.5*(SCLX+SCLY)
                 SCLY=SCLX
              end if
              SCLZ = SCLZ+DQP*(VOL*DPE(4)+2.*TKE/FNST-SCLZ*SC)*HTSTEP
              SCL = SCLX+SCLY+SCLZ   ! trace(P)
            else
              SCL = SCL + DQP*(3.*VOL*DPE(1)+6.*TKE/FNST-SCL*SC)*HTSTEP
              SCLX=SCL
              SCLY=SCL
              SCLZ=SCL
C    Control fluctuations
	      if(abs(SCL*HTSTEP).ge.0.1)then
	        write(*,*)' too strong fluctuations in NPT algorithm'
                if(NUMTASK.gt.1)write(*,*)' at node ',TASKID
	        write(*,*)' scaling factor ',SCL*HTSTEP
   	        if(IPRINT.ge.7)then
	          write(*,*)'------------------->'
	          write(*,*)SCL,TRYCK,TRYCKM
	          write(*,*)VIR1*UNITP/(3.*VOL),VIR2*UNITP/(3.*VOL),
     + PELS1*UNITP/(3.*VOL),PELS2*UNITP/(3.*VOL),VIRB*UNITP/(3.*VOL)
	          write(*,*)'<-------------------'
	        end if
                ISKK=ISKK+1
                if(ISKK.gt.100)then 
	write(*,*)' !!! repeated failure of NPT-algorithm' 
	write(*,*)' !!! restart program with constant volume'
	write(*,*)' !!! or increase thermostat parameter for pressure'
                  call FINAL
                end if
              end if    ! if(LSEP   
            end if
*   3.4.2.3 Calculate scaling coefficients  
            SCVX = 1.-((1.-3./FNST)*SCL+SC)*HTSTEP
            SCVY = SCVX
            SCVZ = SCVX
            DRCX = TSTEP*SCLX
            DRCY = TSTEP*SCLY
            DRCZ = TSTEP*SCLZ
            DRC  = (DRCX+DRCY+DRCZ)/3.
            DRCV = DRC*HTSTEP           !  These are second order corrections 
            DRCX = 1.+DRCX + 0.5*DRCX**2 !  and (may be?) ommitted
            DRCY = 1.+DRCY + 0.5*DRCY**2 !  and (may be?) ommitted
            DRCZ = 1.+DRCZ + 0.5*DRCZ**2 !  and (may be?) ommitted
*
*   3.4.2.4 Correct coordinates and velocities - flexible molecules
            do I=NAB(TASKID),NAE(TASKID)
	        ITYP = ITYPE(I)
	        if(LMOVE(ITYP))then 
                  VX(I) = VX(I)*SCVX
                  VY(I) = VY(I)*SCVY
                  VZ(I) = VZ(I)*SCVZ
                  SX(I) = SX(I)*DRCX + VX(I)*DRCV*MASSDI(I)
                  SY(I) = SY(I)*DRCY + VY(I)*DRCV*MASSDI(I)
                  SZ(I) = SZ(I)*DRCZ + VZ(I)*DRCV*MASSDI(I)
                end if
            end do
*   3.4.2.5  Scale simulation box
            call RECLEN(DRCX,DRCY,DRCZ)
         else      !  .not.LNPT
*  3.4.3 Velocities corrections in the NVT ensemble 
           call GETEMP
           if(LSCTD)then
             do I=NAB(TASKID),NAE(TASKID)
	       ITYP = ITYPE(I)
               SCV = 1.-SCM(ITYP)*HTSTEP
	       if(LMOVE(ITYP))then
                 VX(I) = VX(I)*SCV
                 VY(I) = VY(I)*SCV
                 VZ(I) = VZ(I)*SCV
               end if
             end do
           else
             SCV = 1.-SC*HTSTEP
             do I=NAB(TASKID),NAE(TASKID)
	       ITYP = ITYPE(I)
	       if(LMOVE(ITYP))then
                 VX(I) = VX(I)*SCV
                 VY(I) = VY(I)*SCV
                 VZ(I) = VZ(I)*SCV
               end if
             end do
           end if   ! if(LSCTD)
         end if     ! if(LNPT)
*   3.4.4  Absorbed energy
         EABS = EABS - (SC+SCL*(1.-3./FNST))*TKE*TSTEP
       end if       ! if(LNVT)
*
*   3.5  Second integration of Nose-Hoover equations (t+1/2->t+1)
*   -------------------------------------------------------------
      else if(IAL.eq.2)then
        if(LNVT)then
*   3.5.1 Correct velocities due to thermostat
          if(LSCTD)then
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
              SCV = 1.-SCM(ITYP)*HTSTEP
	      if(LMOVE(ITYP))then
                VX(I) = VX(I)*SCV
                VY(I) = VY(I)*SCV
                VZ(I) = VZ(I)*SCV
              end if
            end do
          else
            SCV = 1.-SC*HTSTEP
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
	      if(LMOVE(ITYP))then 
                VX(I) = VX(I)*SCV
                VY(I) = VY(I)*SCV
                VZ(I) = VZ(I)*SCV
              end if
            end do
          end if
*   3.5.2 Correct velocities due to barostat
          if(LNPT)then
            SCVX = 1.-(1.-3./FNOP)*SCL*HTSTEP
            SCVY = SCVX
            SCVZ = SCVX
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
	      if(LMOVE(ITYP))then 
                VX(I) = VX(I)*SCVX
                VY(I) = VY(I)*SCVY
                VZ(I) = VZ(I)*SCVZ
              end if
            end do
*   3.5.3  Recalculate temperature
            call GETEMP 
*   3.5.4  Correct ita
            if(LSEP)then
              SCLX = SCLX+DQP*(VOL*DPE(2)+2.*TKE/FNST-SCLX*SC)*HTSTEP
              SCLY = SCLY+DQP*(VOL*DPE(3)+2.*TKE/FNST-SCLY*SC)*HTSTEP
              if(LHEX)then
                 SCLX=0.5*(SCLX+SCLY)
                 SCLY=SCLX
              end if
              SCLZ = SCLZ+DQP*(VOL*DPE(4)+2.*TKE/FNST-SCLZ*SC)*HTSTEP
              SCL = SCLX+SCLY+SCLZ   ! trace(P)
            else
              SCL = SCL + DQP*(3.*VOL*DPE(1)+6.*TKE/FNST-SCL*SC)*HTSTEP
              SCLX=SCL
              SCLY=SCL
              SCLZ=SCL
C    Control fluctuations
	      if(abs(SCL*HTSTEP).ge.0.1)then
	        write(*,*)' too strong fluctuations in NPT algorithm'
                if(NUMTASK.gt.1)write(*,*)' at node ',TASKID
	        write(*,*)' scaling factor ',SCL*HTSTEP
   	        if(IPRINT.ge.7)then
	          write(*,*)'------------------->'
	          write(*,*)SCL,TRYCK,TRYCKM
	          write(*,*)VIR1*UNITP/(3.*VOL),VIR2*UNITP/(3.*VOL),
     + PELS1*UNITP/(3.*VOL),PELS2*UNITP/(3.*VOL),VIRB*UNITP/(3.*VOL)
	          write(*,*)'<-------------------'
	        end if
                ISKK=ISKK+1
                if(ISKK.gt.100)then 
	write(*,*)' !!! repeated failure of NPT-algorithm' 
	write(*,*)' !!! restart program with constant volume'
	write(*,*)' !!! or increase thermostat parameter for pressure'
                  call FINAL
                end if
              end if    ! if(LSEP   
            end if
*   3.5.5  Correct ksi due to barostat
            SC = SC + (SCL**2*RTP-DQT/FNST)*HTSTEP
            if(LSCTD)then
              do I=1,NTYPES
                SCM(I) = SCM(I) + (SCL**2*RTP-DQT/FNST)*HTSTEP
              end do
            end if
          else      !  .not.LNPT
*   recalculate temperature
            call GETEMP
          end if    ! if(LNPT
*   3.5.6 Correction ksi due to thermostat
          if(LSCTD)then
            SCA = 0.
            do I=1,NTYPES
              DTE(I)   = TEMPR(I)/TRTEMP-1.d0          
              SCM(I) = SCM(I) + DTE(I)*DQT*HTSTEP
              SCA = SCA + SCM(I)*EKIN(I) 
            end do
            SC = SCA/TKE
          else 
            DKE  = TEMP/TRTEMP-1.d0
	    SC = SC + DKE*DQT*HTSTEP
          end if
*   3.5.7  Absorbed energy
         EABS = EABS - (SC+SCL*(1.-3./FNST))*TKE*TSTEP
        end if ! if (LNVT
*
*   3.6  Case of single-step SHAKE algorithm
*   ----------------------------------------
      else if(IAL.eq.3)then
*   3.6.1  Correct ksi(t+1) from velocities and temperature at t+1/2
        if(LNVT)then
	  SC = SC + DKE*DQT*TSTEP      !  Ksi(t+1)
C    This is separate thermostat for each species (may be not work)
          if(LSCTD)then
            SCA = 0.
            do I=1,NTYPES
              SCM(I) = SCM(I) + DTE(I)*DQT*TSTEP     !  DQT -thermostat mass
              SCA = SCA + SCM(I)*EKIN(I) 
            end do
            SC = SCA/TKE
          end if
          if(LNPT)then
*   3.6.2 Barostat corrections (NPT)
*   -------------------------------
*   3.6.2.1 Correction ita  t-1/2  -> t+1/2 
C         TKE is kinetic energy
            if(LSEP)then
              SCLX = SCLX+DQP*(VOL*DPE(2)+2.*TKE/FNST-SCLX*SC)*TSTEP
              SCLY = SCLY+DQP*(VOL*DPE(3)+2.*TKE/FNST-SCLY*SC)*TSTEP
              if(LHEX)then
                 SCLX=0.5*(SCLX+SCLY)
                 SCLY=SCLX
              end if
              SCLZ = SCLZ+DQP*(VOL*DPE(4)+2.*TKE/FNST-SCLZ*SC)*TSTEP
              SCL = SCLX+SCLY+SCLZ   ! trace(P)
            else
              SCL = SCL + DQP*(3.*VOL*DPE(1)+6.*TKE/FNST-SCL*SC)*TSTEP
              SCLX=SCL
              SCLY=SCL
              SCLZ=SCL
*   3.6.2.2    Control fluctuations
	      if(abs(SCL*TSTEP).ge.0.1)then
	        write(*,*)' too strong fluctuations in NPT algorithm'
                if(NUMTASK.gt.1)write(*,*)' at node ',TASKID
	        write(*,*)' scaling factor ',SCL*TSTEP
   	        if(IPRINT.ge.7)then
	          write(*,*)'------------------->'
	          write(*,*)SCL,TRYCK,TRYCKM
	          write(*,*)VIR1*UNITP/(3.*VOL),VIR2*UNITP/(3.*VOL),
     + PELS1*UNITP/(3.*VOL),PELS2*UNITP/(3.*VOL),VIRB*UNITP/(3.*VOL)
	          write(*,*)'<-------------------'
	        end if
                ISKK=ISKK+1
                if(ISKK.gt.100)then 
	write(*,*)' !!! repeated failure of NPT-algorithm' 
	write(*,*)' !!! restart program with constant volume'
	write(*,*)' !!! or increase thermostat parameter for pressure'
                  call FINAL
                end if
              end if    ! if(LSEP   
            end if
*  this is for correction of velocities
            SCV = ((1.-3./FNST)*SCL)*TSTEP
*   3.6.2.3 Calculate scaling coefficients  
            DRCX = TSTEP*SCLX       ! eta*dt
            DRCY = TSTEP*SCLY
            DRCZ = TSTEP*SCLZ
*   3.6.2.4 Correct coordinates and velosities 
*   3.6.2.4.1   Calculate centre of mass and molecular momenta
	    IMB = NNUM(NAB(TASKID))
	    IME = NNUM(NAE(TASKID))
C    Cycle over molecules
	    do IMOL = IMB,IME
              ITYP = ITM(IMOL)
              SUMM        = SUMMAS (ITYP)
	      if(LMOVE(ITYP))then 
		NSBEG = ISADR(ITYP)+1
		NSEND = ISADR(ITYP+1)
		XC       = 0.D0
      	        YC       = 0.D0
                ZC       = 0.D0
	        PXC	= 0.
                PYC	= 0.
                PZC	= 0.
C   Cycle over atoms within the molecule                  
		IBEG=ISADDR(ITYP)+1+NSITS(ITYP)*(IMOL-IADDR(ITYP)-1)
		IEND=ISADDR(ITYP)+NSITS(ITYP)*(IMOL-IADDR(ITYP))
                DO I       = IBEG,IEND
                  IS	     = NSITE(I)
                  XC      = XC+MASS(IS)*SX(I)
                  YC      = YC+MASS(IS)*SY(I)
                  ZC      = ZC+MASS(IS)*SZ(I)
C   Molecular momenta
                  PXC     = PXC+VX(I)
                  PYC     = PYC+VY(I)
                  PZC     = PZC+VZ(I)
                END DO! OF N
                XC       = XC/SUMM    
      	        YC       = YC/SUMM
      	        ZC       = ZC/SUMM
*   3.6.2.4.2  Correction to COM momenta
                DVX	= PXC*SCV
                DVY	= PYC*SCV
                DVZ	= PZC*SCV 
*   3.6.2.4.3  COM displacements
                DXC = XC*DRCX
                DYC = YC*DRCY
                DZC = ZC*DRCZ
*   3.6.2.4.4  Correct atom velocities
C    DVX - molecular momentum correction
C    DVX/SUMM - COM velocity correction
C    MASS(IS)*DVX/SUMM - atom momentum correction
      	        DO I       = IBEG,IEND
      	          IS	= NSITE(I)
                  FAC	= MASS(IS)/SUMM
                  VX(I)	= VX(I)-DVX*FAC
                  VY(I)	= VY(I)-DVY*FAC
                  VZ(I)	= VZ(I)-DVZ*FAC
*   3.6.2.4.5  Correct atom positions
                  SX(I) = SX(I) + DXC 
                  SY(I) = SY(I) + DYC 
                  SZ(I) = SZ(I) + DZC 
                end do
	      end if
            end do ! of ITYP
*   3.6.2.5  Scale simulation box
            DRCX=DRCX+1.
            DRCY=DRCY+1.
            DRCZ=DRCZ+1.
            call RECLEN(DRCX,DRCY,DRCZ)
*   3,6,2,6 Correction ksi due to barostat
C
C   Additional correction to ksi:
C       d(ksi)=ita**2/W-kT
C       ita/W -> SCL
C     so d(ksi)  ->   SCL**2*(W/Q) - kT/Q 
C                       
            SC = SC + (SCL**2*RTP-DQT/FNST)*TSTEP    !  RTP=DQT/DQP
            if(LSCTD)then    
              SCA = 0.
              do I=1,NTYPES
                SCM(I) = SCM(I) + (SCL**2*RTP-DQT/FNST)*TSTEP
                SCA = SCA + SCM(I)*EKIN(I) 
              end do
              SC = SCA/TKE
            end if
          end if     ! if(LNPT)
*   3.6.4  Absorbed energy
          EABS = EABS - (SC+SCL*(1.-3./FNST))*TKE*TSTEP
        end if       ! if(LNVT)
      end if   ! if(IAL
*
*  3.7 Forcible adjusting of velocities to given temperature
*  ---------------------------------------------------------
      IF(LSCLT.and.(IAL.eq.2.or.IAL.eq.3)) THEN 
        if(LSCTD)then
          DO ITYP     = 1,NTYPES
            IF(DABS(TRTEMP-TEMPR(ITYP)).GT.TDELT.and.LMOVE(ITYP))THEN
              SCV = DSQRT(TRTEMP/TEMPR(ITYP))
	      IBEG=ISADDR(ITYP)+1
	      IEND=ISADDR(ITYP+1) 
	      if(IBEG.le.NAB(TASKID))IBEG=NAB(TASKID)
	      if(IEND.gt.NAE(TASKID))IEND=NAE(TASKID)
	      do I=IBEG,IEND
                VX(I)       = SCV*VX(I)
                VY(I)       = SCV*VY(I)
                VZ(I)       = SCV*VZ(I) 
	      end do  
              if(IPRINT.ge.5)write(*,*)
     +      'velocities scaled by ',SCV,' for type',ITYP
              NRTSC(ITYP) =  NRTSC(ITYP)+1
	    end if
	  end do
        else
          IF(DABS(TRTEMP-TEMP).GT.TDELT)THEN
            SCV = dsqrt(TRTEMP/TEMP)
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
	      if(LMOVE(ITYP))then 
                VX(I) = VX(I)*SCV
                VY(I) = VY(I)*SCV
                VZ(I) = VZ(I)*SCV
              end if
            end do
            if(IPRINT.ge.5)write(*,*)
     +      'velocities scaled by ',SCV
            NRTSC(ITYP) =  NRTSC(ITYP)+1
          end if
        end if
        call GETEMP
      END IF!(LSC...
*
      RETURN
      END    
*
*================== PI_AVER ======================================
*
      subroutine PI_AVER
      include "common_inc.f"
      include "mpif.h"
      real*8 BUFF(NBUFF),BUF2(NBUFF)
      BUFF(1) = PELS1
      BUFF(2) = PELS2
      BUFF(3) = SPE
      BUFF(4) = PE1
      BUFF(5) = PE2
      BUFF(6) = PEST
      BUFF(7) = PINT
      BUFF(8) = TKE
      BUFF(9) = TROT
      BUFF(10) = TTR
      BUFF(11) = TINT
      BUFF(12) = TEMP
      BUFF(13)=VIR1
      BUFF(14)=VIR2
      BUFF(15)=VIRB
      BUFF(16)=WIRSUM
      BUFF(17) = VIRD
      IBUF=17+MOLINT
      do INT=1,MOLINT
	 BUFF(17+INT)=POTLJ(INT)
	 BUFF(IBUF+INT)=POTES(INT)
      end do
      IBUF=IBUF+MOLINT
      IBUF2=IBUF+NTYPES
      do I=1,NTYPES
	 BUFF(IBUF+I)=PES14(I)
	 BUFF(IBUF2+I)=PSR14(I)
      end do
      IBUF=IBUF2+NTYPES
      if(ISTAV.gt.0)then
	if(NRBB.gt.0)then 
	  do I=1,NRBB
	    BUFF(IBUF+I)=BE(I)
	  end do
          IBUF=IBUF+NRBB
        end if
	if(NRAA.gt.0)then 
	  do I=1,NRAA
	    BUFF(IBUF+I)=AE(I)
	  end do
	  IBUF=IBUF+NRAA
	end if
	if(NRTT.gt.0)then 
	  do I=1,NRTT
	    BUFF(IBUF+I)=TE(I)
	  end do
	  IBUF=IBUF+NRTT
	end if
      end if     ! if(ISTAV
      if(IBUF.gt.NBUFF)then
	write(*,*) 'increase NBUFF to ',IBUF
	call FINAL
      end if       
      call MPI_ALLREDUCE(BUFF,BUF2,IBUF,MPI_DOUBLE_PRECISION,
     +    MPI_SUM,MPI_COMM_WORLD,ierr)
      PELS1 = BUF2(1)/NBEADS
      PELS2 = BUF2(2)/NBEADS
      SPE = BUF2(3)/NBEADS
      PE1 = BUF2(4)/NBEADS
      PE2 = BUF2(5)/NBEADS
      PEST = BUF2(6)/NBEADS
      PINT = BUF2(7)/NBEADS
      TKE = BUF2(8)/NBEADS
      TROR = BUF2(9)/NBEADS
      TTR = BUF2(10)/NBEADS
      TINT = BUF2(11)/NBEADS
      TEMP = BUF2(12)/NBEADS
      VIR1 = BUF2(13)/NBEADS
      VIR2 = BUF2(14)/NBEADS
      VIRB = BUF2(15)/NBEADS
      WIRSUM = BUF2(16)/NBEADS
      VIRD = BUF2(17)/NBEADS
      IBUF=17+MOLINT
      do INT=1,MOLINT
	 POTLJ(INT)=BUF2(INT+17)/NBEADS
	 POTES(INT)=BUF2(IBUF+INT)/NBEADS
      end do
      IBUF=IBUF+MOLINT
      IBUF2=IBUF+NTYPES
      do I=1,NTYPES
	 PES14(I)=BUF2(IBUF+I)/NBEADS
	 PSR14(I)=BUF2(IBUF2+I)/NBEADS
      end do
      IBUF=IBUF2+NTYPES
      if(ISTAV.gt.0)then
	if(NRBB.gt.0)then 
	  do I=1,NRBB
	    BE(I)=BUF2(IBUF+I)/NBEADS
	  end do
          IBUF=IBUF+NRBB
        end if
	if(NRAA.gt.0)then 
	  do I=1,NRAA
	    AE(I)=BUF2(IBUF+I)/NBEADS
	  end do
	  IBUF=IBUF+NRAA
	end if
	if(NRTT.gt.0)then 
	  do I=1,NRTT
	    TE(I)=BUF2(IBUF+I)/NBEADS
	  end do
	  IBUF=IBUF+NRTT
	end if
      end if
      return
      end
