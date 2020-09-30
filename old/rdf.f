*
*================= RDFINP =============================================
*
      SUBROUTINE RDFINP
      include"prcm.h" 
*     
      if(TASKID.eq.MAST)then
        PRINT *
        PRINT *,'----------------------------------------'
        PRINT *,'*** PREPARATION FOR RDF CALCULATIONS ***'
        PRINT *,'----------------------------------------'
      end if
      if(TASKID.eq.MAST)
     +PRINT "('NR OF PAIRS IN RDF CALCULATIONS: ',I5)",MAXRDF
      IF(MAXRDF.GT.MAXGR) STOP '!!! MAXRDF.GT.MAXGR IN RDFINP'
*
      N=MAXRDF
      if(TASKID.eq.MAST)
     +PRINT "('*** NR OF RADIAL DISTRIBUTION FUNCTIONS: ',I4)",N
*
      NRDF       = 0
      DO I       = 1,MAXGR
	  IIN(I)=0
	  IGR(I)=0
	  JGR(I)=0
        DO J       = 1,MAX
          IRDF(J,I)  = 0
        END DO
      END DO
*     
*  Check RDF list and calculate normalising factors
      DO ITYP    = 1,NTYPES
        NSPI       = NSPEC(ITYP)
        ISB        = ISADR(ITYP)+1
        ISE        = ISADR(ITYP +1)
        DO IS      = ISB,ISE 
	    NNR=NGRI(IS)
	    do K=1,NNR
	      JS = IGRI(IS,K)
	      if(IS.le.JS)then
		N    = MGRI(IS,K)
		JTYP = ITS(JS)
              NSPJ = NSPEC(JTYP)
              NMGR(N)    = NM(IS)//NM(JS)
              IGR (N)    = NSPI
              JGR (N)    = NSPJ
	        if(IS.eq.JS.and.NSPI.ne.1)then
	          IIN(N)=IIN(N)+(NSPI-1)*NSPI/2
	        else
	          IIN(N)=IIN(N)+NSPI*NSPJ
	        end if     
	      end if
	    end do
*
        END DO! OF IS
      END DO! OF ITYP
*
      DRDF     = RDFCUT/DFLOAT(MAX)
      DRDFI    =  1.D0/DRDF
      if(TASKID.eq.MAST)then
        PRINT "('*** ONE RDF SLICE IN ANGSTROM:  ',F8.5)",DRDF
        PRINT *
	end if 
*  Check RDF table
	if(IPRINT.ge.7)then   
	  write(*,*)' Calculate RDF:'  
	  do IS=1,NSITES
	    NNR=NGRI(IS) 
          do KS=1,NNR
	      JS=IGRI(IS,KS)            
	      K=MGRI(IS,KS)
	    write(*,'(a,2i5,3x,a4,a1,a4,i6,3i9)')
     +' sites ',IS,JS,NM(IS),'-',NM(JS),K,IGR(K),JGR(K),IIN(K)/2
	    end do                                     
	  end do
	end if
*
      RETURN      
      END                       
*
*================= RDFOUT =============================================
*
      SUBROUTINE RDFOUT(IOUT)
      include"prcm.h"
*
      CHARACTER*12 LABRDF
      DIMENSION GR(MAXS,MAXGR),GRO(MAXS,MAXGR),DF(MAXS)
      integer IRDF1(MAXS,MAXGR)
      DATA JINIT/0/, LABRDF /'rdf:        '/
**      save GRO,FNRNR,JINIT,NRNR
*     
      IF(NRDF.EQ.0.and.(.not.LGRST)) RETURN
CM
CM   collecting RDF
CM    
      if(.not.LPIMD)call ADD_RDF(IRDF,IRDF1,MAXS,MAXGR,MAST)  
      if(TASKID.ne.MAST)RETURN
      if(IOUT.ne.0)then
      PRINT *
      PRINT *,'+++++++++++++++++++++++++++++++++++++++++++++++'
      PRINT *,'***      OUTPUT FROM RDF CALCULATIONS       ***'
      PRINT *,'+++++++++++++++++++++++++++++++++++++++++++++++'
      PRINT *
      PRINT "('*** RDF  COLLECTED  DURING:  ',I6,'   STEPS')",NRDF
      PRINT "('*** NR  OF  RDF TYPES  CALCULATED       ',I6)",MAXRDF
      PRINT "('*** CUT OFF RADIUS DIVIDED TO ',I6,' SLICES')",MAX
      PRINT *
*
      PRINT "(1X,130('='))"
      end if
      FNRDF        = DFLOAT(NRDF)
      if(JINIT.eq.0)then
      NRNR         = 0
      FNRNR        = 0.D0
      DO I         = 1,MAXRDF
      DO J         = 1,MAX
      GRO(J,I)      = 0.D0
      END DO! OF J
      END DO! OF I
      end if
*
      IF(LGRST.and.JINIT.eq.0) THEN
        OPEN (UNIT=20,file=frdf,status='old',FORM='UNFORMATTED',
     +        err=111)
        READ(20,err=111) NRNR
        DO I         = 1,MAXRDF
          READ(20,err=112,end=112) (GRO(J,I),J=1,MAX) 
        END DO
        PRINT "('*** RDF RESTART FOR ',I8,' OLD CONFIGURARTIONS')",NRNR
        FNRNR        = DFLOAT(NRNR)
        close (20)
        JINIT=1
      END IF 
	go to 113
 111  write(*,*)' RDF restart failed. File ',frdf
	write(*,*)' RDF are collected only from this run'
	go to 113
 112	write(*,*)' Invalid RDF restart file ',frdf
	frdf(1:1)='?'
	if(LGDMP)write(*,*)' RDF dump will be in file',frdf 
 113  continue
*
      FNRINV       =  1.D0/(FNRDF+FNRNR)
      if(NRDF.eq.0)then
        FNRDF=1.
      end if
      RDFINC       =  RDFCUT/DFLOAT(MAX)
*
      DO 200 I     = 1,MAXRDF
*
        FNI          = DFLOAT(IGR(I))
        FNJ          = DFLOAT(JGR(I))
	FIJ	     = dfloat(IIN(I)) 
*
      if(IOUT.ne.0)then
        write(*,*)
        write(*,'(a,a8,a,i8)')'#*** THIS PAIR: ',NMGR(I),
     +'   --->  TOTAL NR OF CONFIGURATIONS',NRNR+NRDF
        write(*,*)
        write(*,'(3X,5HR(A) ,5X,3HRDF,3X,7HINT RDF)')
*
      end if
      AVRNRI       =  FNRDF
      DNSNRI       = VOL/FNI
      DNSNRJ       = VOL/FNJ 
      RDFINT       = 0.D0 
      LABRDF(5:12)=NMGR(I)
*
      M            = 0
      DO 100 J     = 1,MAX
        M            = M+1    
        DIST         = RDFINC*(DFLOAT(J)-0.5)
        SHELLV       = 4.D0*PI*RDFINC*(DIST**2+RDFINC**2/12.D0)
        FIRDF        = DFLOAT(IRDF1(J,I))
        RDFINT       = RDFINT+FIRDF/ (AVRNRI*FNI)
        RDFIJ        = VOL*FIRDF/(FIJ*AVRNRI*SHELLV)
*
        GR(J,I)      = (FNRDF*RDFIJ+FNRNR*GRO(J,I))*FNRINV
        RDFIJ        = GR(J,I)
	RDFINT	     = RDFINT+RDFIJ*SHELLV/DNSNRJ
        DF(M)        =  RDFIJ*DIST**2
        FINT         =  4.D0*PI*XINTGR(RDFINC,DF,M)
        FINTI        =  FINT/DNSNRI
        FINTJ        =  FINT/DNSNRJ
*
        IF(RDFINT.NE.0.D0.AND.DIST.LE.RDFCUT.and.IOUT.ne.0)PRINT 
     X"(4(1X,F8.4),10x,a12,I9)",
     +DIST,RDFIJ,FINTI,FINTJ,LABRDF
  100   CONTINUE
  200 CONTINUE
*
      if(IOUT.ne.0)PRINT "(1X,130('='))"       
	NRDFT=NRDF+NRNR
      IF(LGDMP.and.(.not.linput).and.NRDFT.gt.0) THEN
      OPEN (UNIT   = 21,file=frdf,status='unknown',FORM='UNFORMATTED')
      WRITE(21) NRDFT
      DO I         = 1,MAXRDF
      WRITE(21) (GR(J,I),J=1,MAX)
      END DO! OF I
      PRINT "('*** RDF DUMP FOR ',I8,' CONFIGURARTIONS')",NRNR+NRDF
      close (21)
      END IF
*
      RETURN
      END
*
*==================== ADD_RDF ======================================
*                                                                   
	subroutine ADD_RDF(IRDF,IRDF1,MAXS,MAXGR,MAST)
 	integer IRDF(MAXS,MAXGR),IRDF1(MAXS,MAXGR) 
	do I=1,MAXS
	  do J=1,MAXGR
	    IRDF1(I,J)=IRDF(I,J)
	  end do
	end do
	return
	end
*
*
*======================== TRANSRDF ========================
*
 	subroutine TRANSRDF(NS,MXGRS,MAST,ID,NGRI,IGRI,MGRI,IAUX)
	integer NGRI(*),IGRI(*),MGRI(*),IAUX(*)
	return
	end
*
*==================== ADD_RDF ======================================
*                                                                   
	subroutine ADD_RDF(IRDF,IRDF1,MAXS,MAXGR,MAST)
	include "mpif.h"
	integer IRDF(MAXS,MAXGR),IRDF1(MAXS,MAXGR) 
	LEN=MAXS*MAXGR
      call MPI_REDUCE(IRDF,IRDF1,LEN,MPI_INTEGER,MPI_SUM,MAST,
     $     MPI_COMM_WORLD,ierr)                 
	return
	end
