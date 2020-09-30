C========================================
C
C   MDynaMix v.5.2
C
C   Part 10
C
C   File service.f
C  
C   This part contains different auxiliary utilities for MD simulations:
C
C   1. GETEMP - calculate temperagture
C   2. COMVEL - initialize velocities
C   3. SCALEV - scale velocities by given factor
C   4. SCLFRC - cut large forces   
C   5. CHTERM - change temperature or/and density
C   6. ZEROFS - set up zero-th on force arrays
C   7. GETCOM - calculate centre of mass 
C   8. RECLEN - recalculate sizes
C   9. DIPOLE - calculate dipole moment - depreciated (moved to unused.f)
C   10. GETKIN - inertia tenzor and related properties
C   11. GETROT - rotation matrix to principal coord. system
C   12. ROTMAT - recalculate coordinates in principal coord. system
C   13. ZEROAV - put zeros in some arrays
C   14. TESTCONF - check configuration for non-overlapping 
C   15. PBC    - periodic boundary conditions
C   16. GROUP -  try to put the whole molecule in the same periodic cell 
C   17. CHECKMOM - check and eventually remove motion of center-of-mass 
*
*   1. Calculate temperature
*   ------------------------
C   This subroutine calculates total temperature and "temperatures" of
C   each molecular specii.
C=============== GETEMP ==============================================
C
      SUBROUTINE GETEMP
      include "common_inc.f"       
      TKE         = 0.D0
      TEMP        = 0.D0 
      VIRKX = 0.
      VIRKY = 0.
      VIRKZ = 0.
*
*   1.1 Computing temperature
*   ------------------------
C   This is applied if each node knows velocities of only own atoms
*   1.1.1  Calculate kinetic energy of own atoms for this node    
        VSQ         = 0.D0
	do ITYP=1,NTYPES
	  TEMPR(ITYP)=0.
	end do
        DO I     = NAB(TASKID),NAE(TASKID)
	  ITYP=ITYPE(I)
	  if(LMOVE(ITYP))then
*   kinetic contribution to virials over directions (OBS! X,Y,Z virials multiplied by 2)
            VVX = MASSDI(I)*VX(I)**2
            VVY = MASSDI(I)*VY(I)**2
            VVZ = MASSDI(I)*VZ(I)**2
            VSC= VVX+VVY+VVZ
            VIRKX = VIRKX + VVX
            VIRKY = VIRKY + VVY
            VIRKZ = VIRKZ + VVZ
	    VSQ=VSQ+VSC
            TEMPR(ITYP)  = TEMPR(ITYP)+0.5D0*VSC   ! kin. energy
 	  end if ! if(LMOVE...                       
	end do
        TEMP         = 0.5D0*VSQ       ! here: kinetic energy
*   1.1.2  Calculate total kinetic energy and virial from all nodes
        TEMPR(NTYP1) = TEMP
        TEMPR(NTYP1+1)=VIRKX
        TEMPR(NTYP1+2)=VIRKY
        NTYP4=NTYP1+3
        TEMPR(NTYP4)=VIRKZ
C   This is used in the case seperate thermostat for each species
C   All nodes must know kinetic energies for each molecular type 
        call ALLSUM(TEMPR,EKIN,NTYP4)
        TKE = EKIN(NTYP1)
        VIRKX = EKIN(NTYP1+1)
        VIRKY = EKIN(NTYP1+2)
        VIRKZ = EKIN(NTYP4)
*   1.1.3 Calculate temperature
	do ITYP=1,NTYPES
          TEMPR(ITYP)  = EKIN(ITYP)*TFACT(ITYP) 
	end do
	if(IPRINT.ge.8)then
	  do ITYP=1,NTYPES
	    write(*,*)ITYP,TEMPR(ITYP)
	  end do
	end if
*
C  This is the main result:
      TEMP        = TKE*CONVET
CC      write(*,*)' GETEMP:',TKE*PERMOL,TEMP,TEMPR(1),TEMPR(2)
*
*  1.2 Check temperature
*  ---------------------
C  This stop the execution if temperature deviation is too high
C  the program check and prints out system configuration
C  This peace is not executed if IPRINT > 5.
      if(LCHECK_TEMP.and.TASKID.eq.MAST.and.
     +(TEMP.gt.CHKTF*TRTEMP.or.TEMP.le.TRTEMP/CHKTF)
     +.and..not.LPIMD)then 
	call TESTCONF
	if(TASKID.eq.MAST)then
          write(*,*)' Temperature control stopped the program '
           write(*,*)' Temperature: ',TEMP 
	  write(*,*)' System configuration:'
          call CONFDUMP(6)
          call XMOLDUMP
	  write(*,*)' Temperature: ',TEMP 
	  write(*,*)' PROGRAM STOP. '
          write(*,*)
     & ' To continue, switch off temperature controll option'  
	end if       
	call ABORT
      end if
      RETURN
      END
*   
*  2. Initialize velocities
*=============== COMVEL ========================================
*
      SUBROUTINE COMVEL
      include "common_inc.f"
*
      DO N    = 1,NSTOT
        IS = NSITE(N)
        RTEMP   = DSQRT(200.D0*MASS(IS)/TRTEMP)            
        VX(N)   = RTEMP*GAUSS()
        VY(N)   = RTEMP*GAUSS()
        VZ(N)   = RTEMP*GAUSS()
      END DO! OF N
*
      SUMX    = 0.D0
      SUMY    = 0.D0
      SUMZ    = 0.D0
      DO N    = 1,NSTOT
      SUMX    = SUMX+VX(N)
      SUMY    = SUMY+VY(N)
      SUMZ    = SUMZ+VZ(N)
      END DO! OF N
      SUMX    = SUMX/DFLOAT(NSTOT)
      SUMY    = SUMY/DFLOAT(NSTOT)
      SUMZ    = SUMZ/DFLOAT(NSTOT)
*
      DO N    = 1,NSTOT
      VX(N)   = VX(N)-SUMX
      VY(N)   = VY(N)-SUMY
      VZ(N)   = VZ(N)-SUMZ
      END DO! OF N
*
      RETURN
      END
*
*      3. Scale velocities
*=============== SCALEV ==========================================
*
      SUBROUTINE SCALEV(VX,VY,VZ,FTMP,N)
      IMPLICIT real*8 (A-H,O-Z)
*
      DIMENSION VX(*),VY(*),VZ(*)
      DO I  = 1,N
      VX(I) = VX(I)*FTMP
      VY(I) = VY(I)*FTMP
      VZ(I) = VZ(I)*FTMP
      END DO
      END
*
*=============== SCLFRC ===================================
*
*    4. Cut large forces
*    -------------------
C
C    This subroutine scale down large forces arising e.g. if
C    initial conformation was too bad. Input parameter is
C    NNN, which says which forces should be scaled:
C    1 - HX (fast forces in double time step algorithm) 
C    2 - GX (slow forces in double time step algorithm) 
C    3 - FX (forces in shake algorithm)
C    Output parameter DONE=.true. if scaling of forces have been done
C    otherwise it returns .false.
C
      SUBROUTINE SCLFRC(NNN)
      include "common_inc.f"
*
      LSCFL      = .FALSE.
      if(.not.LSCFT)return
      FMAX=1.
*
      GO TO (1000,2000,3000),NNN
*
      STOP '!!! WRONG OPTION IN SCLFRC !!!'
*
*  4.1 Cut large fast forces in double time step algorithm
*
 1000 CONTINUE
      FCT=FTSCL*FTIMES
      DO I      = NAB(TASKID),NAE(TASKID)
        FSTERM    = HX(I)**2+HY(I)**2+HZ(I)**2
        FSTERM    = DSQRT(FSTERM)*FCT*MASSDI(I)
        IF(FSTERM.GT.1.D0) THEN
          if(FSTERM.gt.1.d20)then
            FMAX = 1.d20
            write(*,*)' Overflow. Force on atom ',I,' is zeroed '
            HX(I)=0.
            HY(I)=0.
            HZ(I)=0.
          else
            if(FSTERM.gt.FMAX)FMAX=FSTERM
            if(IPRINN.ge.6)
     + write(*,*)' Force on atom ',I,' scaled by ',FSTERM
            HX(I)     = HX(I)/FSTERM
            HY(I)     = HY(I)/FSTERM
            HZ(I)     = HZ(I)/FSTERM
          end if
          LSCFL      = .TRUE.
        END IF
      END DO! OF I
      if(LSCFL)then
         write(*,'(a,i10,a,i4,a,e12.4)')
     +'Fast forces have been cut at step ',
     +MSTEP,' node ',TASKID,' by ',FMAX
      end if
      RETURN
*
*  4.2 Cut large slow forces in double time step algorithm
*
 2000 CONTINUE
      DO I      = NAB(TASKID),NAE(TASKID)
        FSTERM    = GX(I)**2+GY(I)**2+GZ(I)**2
        FSTERM    = DSQRT(FSTERM)*FTSCL*MASSDI(I)
        IF(FSTERM.GT.1.D0) THEN
          if(FSTERM.gt.FMAX)FMAX=FSTERM
          if(FSTERM.gt.1.d20)then
            FMAX = 1.d20
            write(*,*)'Overflow.  Force on atom ',I,' is zeroed '
            GX(I)=0.
            GY(I)=0.
            GZ(I)=0.
          else
            GX(I)     = GX(I)/FSTERM
            GY(I)     = GY(I)/FSTERM
            GZ(I)     = GZ(I)/FSTERM
          if(IPRINN.ge.6)write(*,*)' Force on atom ',I,' scaled by ',
     +    FSTERM
          end if
          LSCFL      = .TRUE.
        END IF
      END DO! OF I
      if(LSCFL)then
         write(*,'(a,i10,a,i4,a,e12.4)')
     +'Slow forces have been cut at step ',
     +MSTEP,' node ',TASKID,' by ',FMAX
      end if
      RETURN
*
*  4.3 Cut large forces in constraint dynamics algorithm
*
 3000 CONTINUE
      DO I      = NAB(TASKID),NAE(TASKID)
        FSTERM    = FX(I)**2+FY(I)**2+FZ(I)**2
        FSTERM    = DSQRT(FSTERM)*FTSCL*MASSDI(I)
        IF(FSTERM.GT.1.D0) THEN
          if(FSTERM.gt.FMAX)FMAX=FSTERM
          if(FSTERM.gt.1.d20)then
            FMAX = 1.d20
            write(*,*)'Overflow. Force on atom ',I,' is zeroed '
            FX(I)=0.
            FY(I)=0.
            FZ(I)=0.
          else
            FX(I)     = FX(I)/FSTERM
            FY(I)     = FY(I)/FSTERM
            FZ(I)     = FZ(I)/FSTERM
          if(IPRINN.ge.6)write(*,*)' Force on atom ',I,' scaled by ',
     +    FSTERM
          end if
          LSCFL      = .TRUE.
        END IF
      END DO! OF I
      if(LSCFL)then
         write(*,'(a,i10,a,i4,a,e12.4)')
     +' Forces have been cut at step ',
     +MSTEP,' node ',TASKID,' by ',FMAX
      end if
*
      RETURN
      END
*
*  5. Change temperature/density and scale velocities/coordinates 
*=================== CHTERM ===================================
* 
*                                                              
      subroutine CHTERM(BOXLN,BOYLN,BOZLN,IND)
      include "common_inc.f"
*  Check excess momenta
      if (LCMOM.and.IND.eq.0) call CHKMOM
*  Set new temperature
	if(LCHT)then
	  SCT=sqrt(TRTEMP/TEMP)
	  do I=1,NSTOT
	    VX(I)=VX(I)*SCT
            VY(I)=VY(I)*SCT
	    VZ(I)=VZ(I)*SCT
          end do
          SC=0.
          do ITYP=1,NTYPES
            SCM(ITYP)=0.
          end do
	  call GETEMP
	  if(TASKID.eq.MAST)write(*,*)' New temperature set to',TEMP
	end if
*  Set new volume/density
	if(LCHP.or.(LVSLOW.and.IND.eq.1))then
	  IF(RHO.gt.1.d-20) THEN
          IF(BOXLN.le.1.d-20.or.BOYLN.le.1.d-20.or.BOZLN.le.1.d-20)THEN
	      VOL = TOTMAS*1.d27/RHO 
	      VOLO = BOXL*BOYL*BOZL
	      if(LOCT)VOLO=0.5*VOLO
	      SCX = (VOL/VOLO)**(1.0/3.0)
	      SCY=SCX
	      SCZ=SCX
              BOXL        = BOXL*SCX
              BOYL        = BOYL*SCY
              BOZL        = BOZL*SCZ
	      if(TASKID.eq.MAST)
     +      write(*,'(2(a,f10.4))')' New density: ',RHO,' BOXL=',BOXL
          else  ! (BOXLN
	      SCX=BOXLN/BOXL                
	      SCY=BOYLN/BOYL                
	      SCZ=BOZLN/BOZL                
	      BOXL=BOXLN
	      BOYL=BOYLN
	      BOZL=BOZLN
	      VOL=BOXL*BOYL*BOZL   
	      if(LOCT)VOL=0.5*VOL
	      DENS	= TOTMAS*1.d-3/(VOL*UNITL**3)
	      if(IND.eq.0.and.TASKID.eq.MAST)
     +      write(*,'(2(a,3f10.4))')'New BOX sizes ',BOXL,BOYL,BOZL,
     +      ' Density',DENS
	      if(IND.eq.1.and.TASKID.eq.MAST.and.IPRINT.ge.7)
     +      write(*,'(2(a,3f10.4))')'New BOX sizes ',BOXL,BOYL,BOZL,
     +      ' Density',DENS
          ENDIF   ! BOXLN
        ELSE      ! RHO
*Vacuum simulation - arbitrary box size:
          CL        = 1000.D0
	    BOXL=CL
	    BOYL=CL
	    BOZL=CL
	    SCX=1.
	    SCY=1.
	    SCZ=1.
          if(TASKID.eq.MAST)PRINT "(/,'*** VACUUM SIMULATION ***',/)"
        END IF    ! RHO
*  Scaling coordinates
        if(LSHEJK)then
	  N=0
	  DO ITYP  = 1,NTYPES                    ! over types
            ISB      = ISADR (ITYP)+1
            ISE      = ISADR (ITYP +1)
            JBE      = IADDR (ITYP)+1
            JEN      = IADDR (ITYP +1)
            DO J     = JBE,JEN                  ! over molecules
	      DX=(SCX-1.d0)*X(J)
	      DY=(SCY-1.d0)*Y(J)
	      DZ=(SCZ-1.d0)*Z(J)
	      X(J)=X(J)+DX
              Y(J)=Y(J)+DY
              Z(J)=Z(J)+DZ
              DO IS    = ISB,ISE                ! over atoms
                N        = N+1
                SX(N)=SX(N)+DX
                SY(N)=SY(N)+DY
                SZ(N)=SZ(N)+DZ
              END DO! OF IS
            end do ! of J
	  end do ! of ITYP
        else
          do I = 1,NSTOT
            SX(I)=SX(I)*SCX
            SY(I)=SY(I)*SCY
            SZ(I)=SZ(I)*SCZ
          end do
          call GETCOM
        end if
	if(IND.eq.0.and.TASKID.eq.MAST)
     + write(*,*)' Scaling distances by ',SCX,SCY,SCZ
	end if  ! LCHP   
	call RECLEN(ODIN,ODIN,ODIN)
	return
	end
*
*=============== ZEROFS =============================================
*  
*  6.  Set some arrays and variables to 0
*
      SUBROUTINE ZEROFS(NNN)
      include "common_inc.f"

*Zero forces:
*
      GO TO (1000,2000),NNN
*
 1000 CONTINUE
      PELS2=0.
      PINT	= 0.
      PE2        = 0.D0
      VIRB=0.
      VIR2=0.
      WIRSS=0.
      VIRFX=0.
      VIRFY=0.
      VIRFZ=0.
      DO I       = 1,MOLINT
	POTE1(I)=0.
	POTL1(I)=0.
      END DO
      do I=1,NTYPES
	PES141(I)=0.d0
	PSR141(I)=0.d0
      end do
      DO I       = 1,NSTOT
        HX(I)      = 0.D0
        HY(I)      = 0.D0
        HZ(I)      = 0.D0
      END DO
      RETURN
*
 2000 CONTINUE
*
      DO I       = 1,MOLINT
	POTE1(I)=0.
	POTL1(I)=0.
      END DO
      do I=1,NTYPES
	PES141(I)=0.d0
	PSR141(I)=0.d0
      end do
      WIRS=0.
C   This collects contributions to virial from intermolecular forces
      VIR1=0.
      PE1=0.
      PELS1=0.
      VIRX=0.
      VIRY=0.
      VIRZ=0.
      DO I       =1,NSTOT
        GX(I)      = 0.D0
        GY(I)      = 0.D0
        GZ(I)      = 0.D0
      END DO
      RETURN
      END
*
*=============== GETCOM ============================================
*                                      
*  7.  Apply periodic boundary conditions for molecular centres of mass
*  --------------------------------------------------------------------
C     The periodic boundary conditions are applied to the molecular
C     centers of mass. The whole molecule is held from one side
C     of the cell (it is not disrupted by PBC)
C
      SUBROUTINE GETCOM
      include "common_inc.f"
      dimension XPC(NPART),YPC(NPART),ZPC(NPART)
*
*  7.1 Calculate molecular centres of mass
*
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES                ! over types
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        SUMM        = SUMMAS (ITYP)
        DO J        = 1,NSPEC(ITYP)         ! over molecules
          I           = I+1
          X (I)       = 0.D0
          Y (I)       = 0.D0
          Z (I)       = 0.D0 
	  PX (I)       = 0.D0
          PY (I)       = 0.D0
          PZ (I)       = 0.D0
C    Calculate C.O.M. vectors for each molecule:
          DO IS       = NSBEG,NSEND          ! over sites
            N           = N+1
            X (I)       = X(I)+MASS(IS)*SX(N)
            Y (I)       = Y(I)+MASS(IS)*SY(N)
            Z (I)       = Z(I)+MASS(IS)*SZ(N)
            PX (I)       = PX(I)+VX(N)
            PY (I)       = PY(I)+VY(N)
            PZ (I)       = PZ(I)+VZ(N)
CD            if(NNUM(N).ne.I)
CD     +write(*,*)' wrong atom/mol in GETCOM -1: ',N,NNUM(N),I,IS,TASKID
          END DO! OF IS
          X (I)       = X(I)/SUMM
          Y (I)       = Y(I)/SUMM
          Z (I)       = Z(I)/SUMM
        END DO! OF I
      END DO! OF ITYP
*
*   7.2 Calculate local atom coordinates relative to molecular COM:
*   --------------------------------------------------------------
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        DO J        = 1,NSPEC(ITYP)
          I           = I+1
          DO IS       = NSBEG,NSEND
            N           = N+1
Calculate short site vectors
            WX(N)       = SX(N)-X(I)
            WY(N)       = SY(N)-Y(I)
            WZ(N)       = SZ(N)-Z(I)
CD            if(NNUM(N).ne.I)
CD     +write(*,*)' wrong atom/mol in GETCOM -2: ',N,NNUM(N),I,IS,TASKID
          END DO! OF IS
        END DO! OF I
      END DO! OF ITYP
*
*   7.3 Apply periodic boundaries on COM points
*   -------------------------------------------
      if(LPIMD)then
        call PI_COM(XPC,YPC,ZPC)           ! COM of centroid 
        do I=1,NOP
          XPC0=XPC(I)
          YPC0=YPC(I)
          ZPC0=ZPC(I)
          call PBC(XPC0,YPC0,ZPC0)
          X(I)=X(I)+XPC0-XPC(I)
          Y(I)=Y(I)+YPC0-YPC(I)
          Z(I)=Z(I)+ZPC0-ZPC(I)
        end do
      else
        DO I        = 1,NOP
          call PBC(X(I),Y(I),Z(I))
        END DO! OF I
      end if  
*
*   7.4  Recalculate new atom coordinates after applying PBC
*   --------------------------------------------------------
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        DO J        = 1,NSPEC(ITYP)
          I           = I+1
          DO IS       = NSBEG,NSEND
            N           = N+1
*Reconstruct long site vectors
            SX(N)       = WX(N)+X(I)
            SY(N)       = WY(N)+Y(I)
            SZ(N)       = WZ(N)+Z(I)
CD            if(NNUM(N).ne.I)
CD     +write(*,*)' wrong atom/mol in GETCOM -3: ',N,NNUM(N),I,IS,TASKID
          END DO! OF IS
        END DO! OF I
      END DO! OF ITYP
*
      RETURN
      END
*
*=============== RECLEN ===============================================
*
C     8. Recalculate sizes
C     --------------------
C     This recalculate all box-size dependent parameters
C     SCX,SCY,SCZ - scaling parameters for each direction
C
      SUBROUTINE RECLEN(SCX,SCY,SCZ)
*
      include "common_inc.f"
*                
      BOXL=SCX*BOXL
      if(LHEX)then
	BOYL=BOXL*cos30  
	BOXY3=0.66666666666666666667d0*BOYL
	BOZL=SCZ*BOZL      
      else if(LOCT)then 
        BOXL=(SCX+SCY+SCZ)*BOXL/(3.*SCX)
	BOYL=BOXL
	BOZL=BOXL
      else
	BOYL=SCY*BOYL
	BOZL=SCZ*BOZL      
      end if
      HBOXL=0.5*BOXL
      HBOYL=0.5*BOYL
      HBOZL=0.5*BOZL 
      RF          = HBOXL
*      RCUTT       = DMIN1(RF,RCUT)
	RCUTT=RCUT
      IF(RCUT.LE.0.D0)then
	  RCUTT = HBOXL
	  if(.not.LHEX.and.RCUTT.gt.HBOYL)RCUTT=HBOYL
	  if(RCUTT.gt.HBOZL)RCUTT=HBOZL
	end if
      RSSQ        = RCUTT**2
      VOL	= BOXL*BOYL*BOZL
	if(ICELL.eq.1)VOL=0.5*VOL
	RHO	= TOTMAS*1.d27/VOL
      return
      end  
*	
*=============== GETKIN ================================================
*
*   10. Calculate inertia tenzor and related properties
*   ---------------------------------------------------
      SUBROUTINE GETKIN
      include "common_inc.f"
*
      real*8 UX(NTOT),UY(NTOT),UZ(NTOT)
      real*8 TOR(9),TORI(9),ROTMX(3,3) 
*                    
      N      = 0
      EKTR	= 0.
      EKROT	= 0.
      TKE	= 0.
      DO ITYP  = 1,NTYPES                    ! over types
	SUMMD	= SUMMAS(ITYP)/(AVSNO*1000.0*UNITM*NBEADS)
        NSS      = NSITS (ITYP)
        NSP      = NSPEC (ITYP)
        ISB      = ISADR (ITYP)+1
        ISE      = ISADR (ITYP +1)
        JBE      = IADDR (ITYP)+1
        JEN      = IADDR (ITYP +1)
        DO J     = JBE,JEN                  ! over molecules
          PX(J)    = 0.D0
          PY(J)    = 0.D0
          PZ(J)    = 0.D0
	  do II=1,9
	    TOR(II)=0.d0
	  end do
          DO IS    = ISB,ISE                ! over atoms
            N      = N+1
*
*	Molecular linear momenta:
*
            PX(J)  = PX(J)+ VX(N)
            PY(J)  = PY(J)+ VY(N)
            PZ(J)  = PZ(J)+ VZ(N)
	    TKE = TKE+0.5*(VX(N)**2+VY(N)**2+VZ(N)**2)*MASSDI(N)
*
*	Molecular angular momenta:
*
	    TOR(1)=TOR(1)+MASSD(IS)*(WY(N)**2+WZ(N)**2)
	    TOR(5)=TOR(5)+MASSD(IS)*(WX(N)**2+WZ(N)**2)
	    TOR(9)=TOR(9)+MASSD(IS)*(WX(N)**2+WY(N)**2)  
	    TOR(2)=TOR(2)-MASSD(IS)*WX(N)*WY(N)
	    TOR(3)=TOR(3)-MASSD(IS)*WX(N)*WZ(N)
	    TOR(6)=TOR(6)-MASSD(IS)*WY(N)*WZ(N)
          END DO! OF IS 
	  TOR(4)=TOR(2)
	  TOR(7)=TOR(3)
	  TOR(8)=TOR(6) 
	  call INV3B3(TOR,TORI,DET,IER)
*  rotation matrix to principal coord.
          if(LMOl1.or.LMOL2)call GETROT
	  EKTR = EKTR+0.5*(PX(J)**2+PY(J)**2+PZ(J)**2)/SUMMD
          RXJ    = 0.D0
          RYJ    = 0.D0
          RZJ    = 0.D0 
	  N	= N-NSS
          DO IS    = ISB,ISE
            N      = N+1
            FAC=MASS(IS)/SUMMAS(ITYP)
* local atomic momenta
            VTX = VX(N)-PX(J)*FAC
            VTY = VY(N)-PY(J)*FAC
            VTZ = VZ(N)-PZ(J)*FAC
* local angular molecular momenta
            RXJ    = RXJ+(VTY*WZ(N)-WY(N)*VTZ)
            RYJ    = RYJ+(VTZ*WX(N)-WZ(N)*VTX)
            RZJ    = RZJ+(VTX*WY(N)-WX(N)*VTY)
          END DO! OF IS    
* rotational kinetic energy
	  if(IER.eq.0)then
	    EROTM=0.5*(TORI(1)*RXJ**2+(TORI(2)+TORI(4))*RXJ*RYJ+
     +                 TORI(5)*RYJ**2+(TORI(3)+TORI(7))*RXJ*RZJ+
     +                 TORI(9)*RZJ**2+(TORI(6)+TORI(8))*RYJ*RZJ)
	  else
	    ROTM=0.
            PP=RXJ**2+RYJ**2+RZJ**2
            if(PP.lt.1.d-50)go to 100
	      N	= N-NSS
            DO IS    = ISB,ISE
              N        = N+1
              SS=WX(N)**2+WY(N)**2+WZ(N)**2
              SP=WX(N)*RXJ+WY(N)*RYJ+WZ(N)*RZJ
              ROTM=ROTM+MASSD(IS)*(SS*PP-SP**2)
            end do
            if(ROTM.ge.1.d-50)then
	        EROTM=0.5d0*PP**2/ROTM
	    else
	        EROTM=0.d0
	    end if
	  end if
	  EKROT=EKROT+EROTM
          TROTM        = EROTM*ENERF*TEMPF
	  if(IPRINT.ge.10)
     +  write(*,'(2i5,e14.5,f10.3)')J,IER,EROTM,TROTM
 100	  continue
	end do ! of J
      end do ! of ITYP
      EKINT=TKE-EKTR-EKROT
      if(FNOP.gt.1.1)then
        TTR         = EKTR*ENERF*TEMPF/(FNOP-1.d0)
      else
        TTR         = EKTR*ENERF*TEMPF/FNOP
      end if
      if(FNSTR.gt.0.01)then
          TROT        = 3.*EKROT*ENERF*TEMPF/FNSTR
	else
	  TROT=TEMP
	end if 
      if(FNSTI.gt.0.01)then
	  TINT	= 3.d0*EKINT*ENERF*TEMPF/FNSTI
	else
	  TINT=0.
	end if
*	TTTT= CONVET*TKE
*	write(*,'(a5,4f12.4)')'kien:',TKET*PERMOL,EKTR*PERMOL,
*     +EKROT*PERMOL,EKINT*PERMOL
*	write(*,'(a5,4f12.4)')'temp:',TTTT,TTR,TROT,TINT
	RETURN
      END
*
*=============== GETROT ==========================================
*
*  11. Calculate rotation matrix to principal coord. system
*  --------------------------------------------------------
      SUBROUTINE GETROT
*
      include"common_inc.f"
      DIMENSION RX(NS),RY(NS),RZ(NS),FMASS(NS)
      DIMENSION FMOI(3,3)      
      DIMENSION G(3),H(3)
*
      N        = 0
      DO ITYP  = 1,NTYPES
        NSS      = NSITS (ITYP)
        NSP      = NSPEC (ITYP)
        ISB      = ISADR (ITYP)+1
        ISE      = ISADR (ITYP +1)
        JBE      = IADDR (ITYP)+1
        JEN      = IADDR (ITYP +1)
        DO J     = JBE,JEN
          I        = 0
          DO IS    = ISB,ISE
            I        = I+1
            N        = N+1
            RX   (I) = WX(N)
            RY   (I) = WY(N)
            RZ   (I) = WZ(N)
            FMASS(I) = MASSD(IS)
          END DO! OF IS
*  Calculate inertia tenzor
          CALL GETMOI(RX,RY,RZ,FMASS,FMOI,NSS)    ! util.f
*  Diagonalize and calculate rotation matrix
          CALL HH3BY3(FMOI,G,H)                   ! util.f
*  Order principal axis according to eigen values of the inertia tenzor
          call ORDER3X3(FMOI,G)                   ! util.f
**        GG=1.d45*UNITM*UNITL**2
**        write(*,'(a,I5,a,3f12.4)')
**     & ' Mol ',I,'  I xx,yy,zz :',GG*G(1),GG*G(2),GG*G(3)
*
          DO K       = 1,3
            DO L       = 1,3
              RMX(J,K,L) = FMOI(K,L)
            END DO! OF K
          END DO! OF L
*
        END DO! OF J
*
      END DO! OF ITYP
*
      RETURN
      END
*
*========= ROTMAT =================================================
*
*   12. Recalculate coordinates in principal coord. system
*   ------------------------------------------------------
      SUBROUTINE ROTMAT(XX,YY,ZZ,IBG,N,NBTOP)
*
      include "common_inc.f"
*
      DIMENSION XX(*),YY(*),ZZ(*)
      LOGICAL NBTOP
*
      IBEG  = IBG
      IEND  = IBG+N-1
*
      IF(NBTOP) THEN	! from box to principal:
        DO I  = IBEG,IEND
          AQ    = XX(I)
          BQ    = YY(I)
          CQ    = ZZ(I)
          XX(I) = RMX(I,1,1)*AQ+RMX(I,2,1)*BQ+RMX(I,3,1)*CQ
          YY(I) = RMX(I,1,2)*AQ+RMX(I,2,2)*BQ+RMX(I,3,2)*CQ
          ZZ(I) = RMX(I,1,3)*AQ+RMX(I,2,3)*BQ+RMX(I,3,3)*CQ
        END DO! OF I
      ELSE 		! from principal to box:
        DO I  = IBEG,IEND
          AQ    = XX(I)
          BQ    = YY(I)
          CQ    = ZZ(I)
          XX(I) = RMX(I,1,1)*AQ+RMX(I,1,2)*BQ+RMX(I,1,3)*CQ
          YY(I) = RMX(I,2,1)*AQ+RMX(I,2,2)*BQ+RMX(I,2,3)*CQ
          ZZ(I) = RMX(I,3,1)*AQ+RMX(I,3,2)*BQ+RMX(I,3,3)*CQ
        END DO! OF I
      END IF
*
      RETURN
      END
*
*=============== ZEROAV =============================================
*
*   13.  Put zeroth in accumulation arrays for averages
*   ---------------------------------------------------
      SUBROUTINE ZEROAV
      include "common_inc.f"
*
*Zero averages
*
      DO I       = 1,NRBB
        BB(I)      = 0.D0
        EB(I)      = 0.D0
      END DO
*
      DO I       = 1,NRAA
        AA(I)      = 0.D0
        EA(I)      = 0.D0
      END DO
*
      DO I       = 1,NRTT
        TT(I)      = 0.D0
        ET(I)      = 0.D0
      END DO 
*
      DO I       = 1,MOLINT
        POTES(I)   = 0.D0
        POTLJ(I)   = 0.D0
      END DO
      do I=1,NTYPES
	SELFPE(I)      = 0.D0
        PES14 (I)      = 0.D0
        PSR14 (I)      = 0.D0
      end do
      RETURN
      END

*
*======================== TESTCONF ===================================
*
*   14. Test configuration for non-overlapping
*   ------------------------------------------
	subroutine TESTCONF
*  Check configuration
	include "common_inc.f"
	logical OK,OKB
	OK=.true.     
	OKB=.true.
	if(TASKID.ne.MAST)return
	do I=1,NSTOT-1 
	  IMOL=NNUM(I)
	  ISB =NSITE(I)
          IST=INBT(ISB)
	  ITYP=ITYPE(I)
	  do J=I+1,NSTOT
	    JMOL=NNUM(J)
	    JSB =NSITE(J)
            JST=INBT(JSB)
	    JTYP=ITYPE(J)
	    DX=SX(I)-SX(J)
	    DY=SY(I)-SY(J)
	    DZ=SZ(I)-SZ(J)
	    call PBC (DX,DY,DZ)
	    RR=sqrt(DX**2+DY**2+DZ**2)
	    if(RR.le.0.25*(SIGT(IST)+SIGT(JST)).or.RR.le.1.)then  
	      if(IMOL.eq.JMOL)then
*  Check that these are non bound
	        do K=1,NNBB(I) 
*  bound
		  if(INBB(I,K).eq.J)go to 125
                  if(INBB(I,K).eq.-J)then
	    if(RR.ge.0.2*(SIGT(IST)+SIGT(JST)).and.RR.ge.1.)go to 125  
                  end if
	  	end do
              end if
*  Check that atoms do not belong EE molecules
              if(IEE(ITYP)*IEE(JTYP).eq.0)go to 125 
              if(OK)then
		  OK=.false.
		  write(*,*)'!!! BAD CONFIGURATION'
		  write(*,*)'The following atom pair are too close'
		  write(*,*)
     +' N1  mol1  site1   N2  mol2  site2    RR        Eel        Elj'
		end if
  	      QIJ=COULF*CHARGE(ISB)*CHARGE(JSB)
              A6        =  A6LJ(IST,JST)
              B12       =  B12LJ(IST,JST)
		EEL=0.001*ENERF*QIJ/RR
		ELJ=0.001*ENERF*(A6/RR**6+B12/RR**12)
		write(*,'(2(I5,1x,a6,1x,a4),f10.4,2e12.4)')
     +I,NAME(ITYP),NM(ISB),J,NAME(JTYP),NM(JSB),RR,EEL,ELJ 
	      end if
 125	    continue
	  end do
	end do
*   Bonds
	do ITYP=1,NTYPES
        NBBEG     = IADB(ITYP)
        NBEND     = NBBEG+NRB(ITYP)-1   
	  do IIS=NBBEG,NBEND
	    do IM=1,NSPEC(ITYP)      
	      ISHIFT=ISADDR(ITYP)+(IM-1)*NSITS(ITYP)
	      I=ISHIFT+IB(IIS)
	      J=ISHIFT+JB(IIS)           
C	write(*,*)' check ',ITYP,IIS,IM,IB(IIS),JB(IIS),I,J
            DX      = SX(I)-SX(J)
            DY      = SY(I)-SY(J)
            DZ      = SZ(I)-SZ(J)
	      call PBC(DX,DY,DZ)
            RR      = sqrt(DX*DX+DY*DY+DZ*DZ) 
	      REQ     = RB(IIS)
	      if(RR.le.0.7*REQ.or.RR.gt.1.5*REQ)then
	        if(OKB)then
	          OKB=.false.
	          write(*,*)'!!! BAD BOND LENGTH:'
	          write(*,*)
     +' Type  Mol     site1       site2      RR      Req        Ebond'
	        end if
	        ISB = ISADR(ITYP)+IB(IIS)
	        JSB = ISADR(ITYP)+JB(IIS)
	        FF=FB(IIS)*(RR-REQ)**2/EFACT
                write(*,'(a6,i5,2(2x,a4,i5),2f10.4,e13.5)')
     +        NAME(ITYP),IM,NM(ISB),ISB,NM(JSB),JSB,RR,REQ,FF
	      end if
	    end do
	  end do
	end do
	if(OK.and.OKB)write(*,*)' Configuration OK'
	return
	end
*
*======================= PBC ========================================
* 
*    15. Periodic boundary conditions
*    --------------------------------
      subroutine PBC(XX,YY,ZZ)
      include "common_inc.f"
      if(abs(XX) > HBOXL) XX = XX - sign(BOXL,XX)
      if(abs(ZZ) > HBOZL) ZZ = ZZ - sign(BOZL,ZZ)
      if(LHEX)then
C  hexagonal periodic cell - minimum image up to 1.5*BOXL from the center    
	XY=BOXYC*XX
        if(YY.gt.BOXY3-XY.and.(XX.gt.0..or.YY.gt.2.*BOXY3-XY))then
          YY=YY-BOYL
          XX=XX-HBOXL
	  XY=BOXYC*XX
        end if
        if(YY.lt.-BOXY3-XY.and.(XX.lt.0..or.YY.lt.-2.*BOXY3-XY))then
          YY=YY+BOYL
          XX=XX+HBOXL
	  XY=BOXYC*XX
        end if
        if(YY.gt.BOXY3+XY)then
          YY=YY-BOYL
          XX=XX+HBOXL
	  XY=BOXYC*XX
        end if
        if(YY.lt.-BOXY3+XY)then
          YY=YY+BOYL
          XX=XX-HBOXL
        end if
      else
        if(abs(YY) > HBOYL) YY = YY - sign(BOYL,YY)
      end if
      if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	CORSS=HBOXL*int((4./3.)*(abs(XX)+abs(YY)+abs(ZZ))/BOXL)
	XX=XX-sign(CORSS,XX)
	YY=YY-sign(CORSS,YY)
	ZZ=ZZ-sign(CORSS,ZZ)
      end if
      return
      end 
*
*================== GROUP ======================================
* 
C  16.  Gather all atoms of a molecule in the same place
C
C
      subroutine GROUP
      include "common_inc.f"
      integer IMASK(NS),JMASK(NS)
      if(TASKID.eq.MAST)
     +write(*,*)' Gathering atoms in each molecule...'
      call GETCOM
      DO ITYP     = 1,NTYPES                ! over types
        if(IINIT(ITYP).ne.1)then
*         write(*,*)' type ',ITYP,NRB(ITYP)
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        SUMM        = SUMMAS (ITYP)
	JSFB=IADB(ITYP)
	JSFE=JSFB+NRB(ITYP)-1
        DO JMOL        = 1,NSPEC(ITYP)         ! over molecules
          I0 = ISADDR(ITYP)+(JMOL-1)*NSITS(ITYP)
          ISHFT = I0 -ISADR(ITYP)
          IS = NSBEG
          do I=1,NSITES
            IMASK(I)=0   ! - atom is OK
            JMASK(I)=0   ! - all neighbouring atoms to this OK
          end do
*          write(*,*)' molecule ',JMOL
          IMASK(IS)=1
C---> start gathering the molecule here
C     curret site is IS and current atom is ISP
C     checking bonds
 10       ISP = IS + ISHFT
*          write(*,*)' proceeding ',IS,ISP,JSFB,JSFE
          do IBB = JSFB,JSFE
            I1 = IB(IBB)+ISADR(ITYP)
            J1 = JB(IBB)+ISADR(ITYP)
            if(IS.eq.I1)then
              if(IMASK(J1).eq.0)then
C connecting atom J1 to I1
                JSP = J1 + ISHFT
                DX = SX(JSP)-SX(ISP)
                DY = SY(JSP)-SY(ISP)
                DZ = SZ(JSP)-SZ(ISP)
                call PBC(DX,DY,DZ)
                call PBC(DX,DY,DZ)
                SX(JSP) = SX(ISP) +DX
                SY(JSP) = SY(ISP) +DY
                SZ(JSP) = SZ(ISP) +DZ
                IMASK(J1)=1
*                write(*,*)' connecting ',JSP,' to ',ISP,' mark ',J1
              end if
            end if
            if(IS.eq.J1)then
              if(IMASK(I1).eq.0)then
C connecting atom J1 to I1
                JSP = I1 + ISHFT
                DX = SX(JSP)-SX(ISP)
                DY = SY(JSP)-SY(ISP)
                DZ = SZ(JSP)-SZ(ISP)
                call PBC(DX,DY,DZ)
                call PBC(DX,DY,DZ)
                SX(JSP) = SX(ISP) +DX
                SY(JSP) = SY(ISP) +DY
                SZ(JSP) = SZ(ISP) +DZ
                IMASK(I1)=1
*                write(*,*)' connecting ',JSP,' to ',ISP,' mark ',I1
              end if
            end if
          END DO! OF IBB
          JMASK(IS)=1
* Finding next OK atom
          do J=NSBEG+1,NSEND
            if(JMASK(J).eq.0.and.IMASK(J).eq.1)then
              IS=J
              go to 10
            end if
          end do
        end do
        end if
      end do
      call GETCOM
      return
      end
*
*================== CHKMOM ==========================================
*
*    17   Check and remove excess COM momenta
*    ----------------------------------------
      subroutine CHKMOM
      include "common_inc.f"
      dimension PT(3),SPT(3)
      do I=1,3
        PT(I)=0.
      end do
*  Check excess momenta
      do I=NAB(TASKID),NAE(TASKID)
	PT(1)=PT(1)+VX(I)
	PT(2)=PT(2)+VY(I)
	PT(3)=PT(3)+VZ(I)
      end do	           
      call GETEMP
      call ALLSUM(PT,SPT,3)
      TEMPO=TEMP
      SUMASD=TOTMAS/UNITM
      DVX=SPT(1)/SUMASD
      DVY=SPT(2)/SUMASD
      DVZ=SPT(3)/SUMASD
      if(TASKID.eq.MAST.and.IPRINT.ge.6)
     + write(*,'(a14,3d18.10)')' COM momenta: ',DVX,DVY,DVZ
      do I=NAB(TASKID),NAE(TASKID)
	VX(I)=VX(I)-DVX/MASSDI(I)
	VY(I)=VY(I)-DVY/MASSDI(I)
	VZ(I)=VZ(I)-DVZ/MASSDI(I)
      end do
*  Check excess momenta
      do I=1,3
        PT(I)=0.
      end do
      do I=NAB(TASKID),NAE(TASKID)
	PT(1)=PT(1)+VX(I)
	PT(2)=PT(2)+VY(I)
	PT(3)=PT(3)+VZ(I)
      end do	           
      call GETEMP
      call ALLSUM(PT,SPT,3)
      DVX=SPT(1)/SUMASD
      DVY=SPT(2)/SUMASD
      DVZ=SPT(3)/SUMASD
         if(TASKID.eq.MAST)
     +      write(*,'(a14,3d18.10)')' COM corrected: ',DVX,DVY,DVZ
      call GETEMP
      TERR = dabs(TEMPO-TEMP)/TEMP
      if(TERR.gt.0.01.and.TASKID.eq.MAST)then
	write(*,*)' anomalous excess momenta removed'
	write(*,*)'old temperature ',TEMPO
        write(*,*)'new temperature ',TEMP
      else
        if(TASKID.eq.MAST)
     +  write(*,'(a,d12.5)')'  Total momentum OK.  Correction ',TERR
      end if
      return
      end




