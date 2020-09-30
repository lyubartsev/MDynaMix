C      Part 6c
C      -------      
C         fake MPI subroutines (for sequential execution)
C         Most of the procedures do nothing.
C         Only some of them make accumulation of data
C         See explanations in mpi.f file
	subroutine PARAINI(NUMTASK,TASKID)
	integer NUMTASK,TASKID
	NUMTASK=1
	TASKID=0
	return
	end
*
*================ SINHR ============================================
*                                                                   
	subroutine SINHR(I)
	return
	end
*
*================ SINHV ============================================
*                                                                   
	subroutine SINHV(I)
	return
	end
*
*=================== CM_ADDF ======================================
*
      subroutine CM_ADDF(LAVER)
      include "common_inc.f"         
      logical LAVER
      timen0	= cputime(0.d0)
      if(LAVER)then
	do INT=1,MOLINT
	  POTES(INT)=POTES(INT)+POTE1(INT)
	  POTLJ(INT)=POTLJ(INT)+POTL1(INT)
	end do           
	do I=1,NTYPES
	  PES14(I)=PES14(I)+PES141(I)
	  PSR14(I)=PSR14(I)+PSR141(I)
	end do
        DO M      = 1,NRBB
          BB(M)     = BB(M)+BR(M)
          EB(M)     = EB(M)+BE(M)
        END DO
        DO M      = 1,NRAA
          AA(M)     = AA(M)+AR(M)
          EA(M)     = EA(M)+AE(M)
	  end do
        DO M      = 1,NRTT
          TT(M)     = TT(M)+TR(M)
          ET(M)     = ET(M)+TE(M)
  	END DO           
      end if    
      timen	= timen+cputime(timen0) 
      end
*
*=================== CM_ADDLF ======================================
*
      subroutine CM_ADDLF
      include "common_inc.f" 
      timen0	= cputime(0.d0)
      do INT=1,MOLINT
	POTES(INT)=POTES(INT)+POTE1(INT)
	POTLJ(INT)=POTLJ(INT)+POTL1(INT)
      end do           
      do I=1,NTYPES
	PES14(I)=PES14(I)+PES141(I)
	PSR14(I)=PSR14(I)+PSR141(I) 
	SELFPE(I)=SELFPE(I)+SPEE(I)
      end do
      DO N         = 1,NSTOT
	M	= NNUM(N)
	WIRS	= WIRS-GX(N)*(SX(N)-X(M))-GY(N)*(SY(N)-Y(M))-
     +               GZ(N)*(SZ(N)-Z(M))
      END DO! OF N
      timen	= timen+cputime(timen0) 
      return
      end
*
*================== CM_ADDLA =============================
*
      subroutine CM_ADDLA
      include "common_inc.f" 
      timen0	= cputime(0.d0)
      do INT=1,MOLINT
	POTES(INT)=POTES(INT)+POTE1(INT)
	POTLJ(INT)=POTLJ(INT)+POTL1(INT)
      end do           
      do I=1,NTYPES
	PES14(I)=PES14(I)+PES141(I)
	PSR14(I)=PSR14(I)+PSR141(I)
	SELFPE(I)=SELFPE(I)+SPEE(I)
      end do
      DO M      = 1,NRBB
        BB(M)     = BB(M)+BR(M)
        EB(M)     = EB(M)+BE(M)
      END DO
      DO M      = 1,NRAA
        AA(M)     = AA(M)+AR(M)
        EA(M)     = EA(M)+AE(M)
      end do
      DO M      = 1,NRTT
        TT(M)     = TT(M)+TR(M)
        ET(M)     = ET(M)+TE(M)
      END DO           
      DO N         = 1,NSTOT
	FX(N)   = GX(N)+HX(N)
	FY(N)   = GY(N)+HY(N)
	FZ(N)   = GZ(N)+HZ(N)
	M	= NNUM(N)
	WIRS	= WIRS-GX(N)*(SX(N)-X(M))-GY(N)*(SY(N)-Y(M))-
     +               GZ(N)*(SZ(N)-Z(M))
      END DO! OF N
      timen	= timen+cputime(timen0) 
      return
      end
*
*===================== GATHER =======================================
*                                                                    
	subroutine GATHER(NRNOD)
	integer NRNOD(*)
	return
	end
*
*============================= ALLSUM ==================================
*     
	subroutine ALLSUM(A,B,N)
	real*8 A(*),B(*)
	do I=1,N
	B(I)=A(I)
	end do
	return
	end
*
*============================ SUMMA ===============================
*                                                                  
      subroutine SUMMA(A,B,N,MAST)
	real*8 A(*),B(*)
	do I=1,N
	B(I)=A(I)
	end do
	return
	end  
*
*============================ BCAST ===============================
*                                                                  
      subroutine BCAST(A,N,MAST)
	real*8 A(*)
	return
	end  
*
*============================ BICAST ===============================
*                                                                  
      subroutine BICAST(I,MAST)
	return
	end  
*
*========================== BARRIER ================================
*
      subroutine BARRIER
      return
      end
*
*============================= FINAL ========================
*                                                            
	subroutine FINAL
	close(6)
	stop 
	end
*
*============================= PIMD ========================
*                                                            
	subroutine PIMD
	stop 
	end
*
*============================= PI_COM =======================
*
	subroutine PI_COM(XPC,YPC,ZPC)
	real*8 XPC(*),YPC(*),ZPC(*)
	return
	end
*
*============================= ABORT ========================
*                                                            
	subroutine ABORT
	stop 
	end
