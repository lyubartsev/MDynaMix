C    Package tranal v.5.0
C    Calculations of radial distribution functions
      program RDF
      include "tranal.h"
      character*64 FOUTRDF
      parameter (MAXRDF=64,MAXGR=1024,NAAM=1024)
      integer*8 IREP(MAXRDF),IRD(MAXRDF,MAXGR),JRD(MAXRDF,MAXGR),
     +IRDF(MAXRDF,NAAM),IRDFI(MAXRDF,NAAM)
      logical LCOMRDF1,LCOMRDF2
      dimension TMP(MAXRDF)
      namelist /RDFIN/ FOUTRDF,RDFCUT,NA,NAI,RMI,RMAX,NRDF,ME,
     +  LCOMRDF1,LCOMRDF2
*  defaults
      LCOMRDF1=.false.
      LCOMRDF2=.false.
      RDFCUT=10.
      NA=200
      NAI=100
      RMI=0
      RMAX=5.
      ME=-1
*  read and setup trajectory parameters
      call SETTRAJ
*  read RDF computation parameters
      read(*,RDFIN)
      if(NRDF.gt.MAXRDF)stop 'Increase MAXRDF!!!'
      if(NA.gt.NAAM.or.NAI.gt.NAAM)stop 'Increase NAAM!'
      DRI=RMAX-RMI
*  initialization
      BOXS=0.           ! to accumulate average box sizes
      BOYS=0.
      BOZS=0.
*  read RDFs
      N=0
      if(LCOMRDF1)write(*,*)
     +' For the first sites in the list, molecular COM is used'
      if(LCOMRDF2)write(*,*)
     +' For the second sites in the list, molecular COM is used'
      write(*,*)' Trying to read ',NRDF,' groups for RDF calculations'
      DO K       = 1,NRDF
 	IE=-1
	STR=TAKESTR(5,IE)
	if(str(1:1).eq.'&')then
          READ(str(2:24),*)IRP
	  if(IRP.gt.MAXGR) stop 'Increase MAXGR!!!'
	  IREP(K)=IRP
	  do IR=1,IRP
	    STR=TAKESTR(5,IE)
	    read(STR,*,err=99)I,J  
	    IRD(K,IR)=I
	    JRD(K,IR)=J
	  end do
	else 
	  IREP(K)=1
	  read(STR,*,err=99)I,J  
	  IRD(K,1)=I
	  JRD(K,1)=J
	end if
	write(*,*)' RDF  No ',K,' equiv ',IREP(K)
        if(LCOMRDF1.and.LCOMRDF2)then
	  do I=1,IREP(K)
	    write(*,'(2I6,2(3x,a6))')
     +      IRD(K,I),JRD(K,I),NAME(ITS(IRD(K,I))),NAME(ITS(JRD(K,I)))
	  end do
        else if(LCOMRDF1.and..not.LCOMRDF2)then
	  do I=1,IREP(K)
	    write(*,'(2I6,2(3x,a6))')
     +      IRD(K,I),JRD(K,I),NAME(ITS(IRD(K,I))),NM(JRD(K,I))
	  end do
        else if(.not.LCOMRDF1.and.LCOMRDF2)then
	  do I=1,IREP(K)
	    write(*,'(2I6,2(3x,a6))')
     +      IRD(K,I),JRD(K,I),NM(IRD(K,I)),NAME(ITS(JRD(K,I)))
	  end do
        else
	  do I=1,IREP(K)
	    write(*,'(2I6,2(3x,a4))')
     +      IRD(K,I),JRD(K,I),NM(IRD(K,I)),NM(JRD(K,I))
	  end do
        end if
      END DO
      go to 101
 99   STR=TAKESTR(5,99)
      stop
 101  continue
      do I=1,NA
	do J=1,NRDF
	  IRDF(J,I)=0
	end do
      end do
      do I=1,NAI
	do J=1,NRDF
	  IRDFI(J,I)=0
	end do
      end do
*  Starting trajectory analysis
      IEND=0
      do while(IEND.eq.0)              ! over configurations
 50     call READCONF(IEND)
        if(LCOMRDF1.or.LCOMRDF2)call GETCOM
        if(IEND.ne.0)go to 200
        if(ME.le.0)then
          if(IPRINT.ge.6)write(*,'(a,f12.3,a,3f10.3)')
     +    ' Processing time',FULTIM
        else
          if(IPRINT.ge.6)write(*,'(a,f12.3,a,i5)')
     +    ' Processing time',FULTIM,'  Ens ',MEE
          if(ME.ne.MEE)then
             IAN=IAN-1
             go to 50
          end if   
        end if
        do I=1,NRDF                    !  over RDF groups
	  do K=1,IREP(I)               !  over site pairs of the same group
	    IS=IRD(I,K)
	    JS=JRD(I,K)
	    ITYP=ITS(IS)
	    JTYP=ITS(JS)
	    NIS=NSPEC(ITYP)
	    NJS=NSPEC(JTYP)
	    do IM=1,NIS                ! over first molecules
	      JBEG=1
	      if(IS.eq.JS)JBEG=IM+1
	      ISHIFT=ISADDR(ITYP)+IS-ISADR(ITYP)-NSITS(ITYP)
	      do JM=JBEG,NJS           ! over second molecules
	        JSHIFT=ISADDR(JTYP)+JS-ISADR(JTYP)-NSITS(JTYP)
	        IST=ISHIFT+IM*NSITS(ITYP)
	        JST=JSHIFT+JM*NSITS(JTYP)
	        IMOL=NNUM(IST)
	        JMOL=NNUM(JST)
                if(LCOMRDF1.and.LCOMRDF2)then
                  DX=X(IMOL)-X(JMOL)
                  DY=Y(IMOL)-Y(JMOL)
                  DZ=Z(IMOL)-Z(JMOL)
                else if(LCOMRDF1.and..not.LCOMRDF2)then
	          DX=X(IMOL)-SX(JST)
	          DY=Y(IMOL)-SY(JST)
	          DZ=Z(IMOL)-SZ(JST)
                else if(.not.LCOMRDF1.and.LCOMRDF2)then
                  DX=SX(IST)-X(JMOL)
                  DY=SY(IST)-Y(JMOL)
                  DZ=SZ(IST)-Z(JMOL)
                else
	          DX=SX(IST)-SX(JST)
	          DY=SY(IST)-SY(JST)
	          DZ=SZ(IST)-SZ(JST)
                end if
	        call PBC(DX,DY,DZ)
	        RR=sqrt(DX**2+DY**2+DZ**2)
**		  write(*,*)I,K,IS,JS,IST,JST,IMOL,JMOL,RR
*  This distinguishes sites on the same molecule
	        if(IMOL.ne.JMOL)then
		  NR=NA*RR/RDFCUT+1
		  if(NR.gt.0.and.NR.le.NA)IRDF(I,NR)=IRDF(I,NR)+1
	        else
		  NR=NAI*(RR-RMI)/DRI+1
**		write(*,'(7i6,f8.4)')I,K,IS,JS,IST,JST,NR,RR
		  if(NR.gt.0.and.NR.le.NAI)IRDFI(I,NR)=IRDFI(I,NR)+1
	        end if
	      end do
	    end do
	  end do
        end do
        BOXS=BOXS+BOXL
        BOYS=BOYS+BOYL
        BOZS=BOZS+BOZL
      end do
 200  write(*,*)IAN,' configurations analysed'
*  RDF output
      open(unit=16,file=FOUTRDF,status='unknown')
      write(16,'(a)')'   Radial distribution functions '
      write(16,*)IAN,' configurations analysed'
*   Averagre box sizes and volume
      BOXS=BOXS/IAN
      BOYS=BOYS/IAN
      BOZS=BOZS/IAN
      VOLS = BOXS*BOYS*BOZS
      write(16,'(a,3f12.4)')'  Average box sizes: ',BOXS,BOYS,BOZS
      IFST=NA
      do K=1,NRDF
	CONU=0.
	NPAIR=0
	NPRI=0                   ! for intramolecular. pairs
	do KK=1,IREP(K)
	  IS=IRD(K,KK)
	  JS=JRD(K,KK)
	  ITYP=ITS(IS)
	  JTYP=ITS(JS)
	  NN1 = NSPEC(ITYP)
	  NN2 = NSPEC(JTYP)
	  if(ITYP.eq.JTYP)then
	    if(IS.eq.JS)then
	      NPAIR=NPAIR+NN1*(NN1-1)/2
	    else
	      NPAIR=NPAIR+NN1*(NN1-1)
	    end if
	  else
	    NPAIR=NPAIR+NN1*NN2
	  end if
	  if(ITYP.eq.JTYP.and.IS.ne.JS)NPRI=NPRI+1
	end do
        write(16,*)
	FACT=VOLS/(NPAIR*1.d0*IAN)
        if(LCOMRDF1.and.LCOMRDF2)then
	  write(16,'(5a,I12,a,I5)')' This mol. pair: ',NAME(ITS(IS)),
     +' - ',NAME(ITS(JS)),'   Npairs=',NPAIR,'  Nint=',NPRI
        else if(LCOMRDF1.and..not.LCOMRDF2)then
	  write(16,'(5a,I12,a,I5)')' This  pair: ',NAME(ITS(IS)),
     +' - ',NM(JS),'   Npairs=',NPAIR,'  Nint=',NPRI
        else if(.not.LCOMRDF1.and.LCOMRDF2)then
	  write(16,'(5a,I12,a,I5)')' This pair: ',NM(IS),
     +' - ',NAME(ITS(JS)),'   Npairs=',NPAIR,'  Nint=',NPRI
        else
	  write(16,'(5a,I12,a,I5)')' This pair: ',NM(IS),' - ',NM(JS),
     +  '   Npairs=',NPAIR,'  Nint=',NPRI
        end if
        write(16,'(a)')
     + '  R          RDF         RCN1        RCN2          KBI(nm^3)'
	write(16,*)'-------------------------------------------------'
        CONU=0.
        RKBI=0.
	do I=1,NA
	  RDFV=IRDF(K,I)
	  RR=(I-0.5)*RDFCUT/NA
	  DRR=RDFCUT/NA
   	  SH12=4*PI*DRR*(RR**2+DRR**2/12.d0)
	  CONU=CONU+RDFV*FACT/VOLS
	  CONU1=CONU*NSPEC(ITS(IS))
	  CONU2=CONU*NSPEC(ITS(JS))
          GofR = RDFV*FACT/SH12
          RKBI=RKBI+(GofR-1.)*SH12*0.001
  	  if(CONU.gt.1.d-18)then
            if(LCOMRDF1.and.LCOMRDF2)then
	      write(16,'(f8.3,3f12.5,2x,f12.5,3x,a4,2x,a4,a1,a4,I12)')
     +RR,GofR,CONU1,CONU2,RKBI,'rdf:',NAME(ITS(IS)),'-',
     +NAME(ITS(JS)),IRDF(K,I)
            else if(LCOMRDF1.and..not.LCOMRDF2)then
	      write(16,'(f8.3,3f12.5,2x,f12.5,3x,a4,2x,a4,a1,a4,I12)')
     +RR,GofR,CONU1,CONU2,RKBI,'rdf:',NAME(ITS(IS)),'-',
     +NM(JS),IRDF(K,I)
            else if(.not.LCOMRDF1.and.LCOMRDF2)then
	      write(16,'(f8.3,3f12.5,2x,f12.5,3x,a4,2x,a4,a1,a4,I12)')
     +RR,GofR,CONU1,CONU2,RKBI,'rdf:',NM(IS),'-',
     +NAME(ITS(JS)),IRDF(K,I)
            else
	      write(16,'(f8.3,3f12.5,2x,f12.5,3x,a4,2x,a4,a1,a4,I12)')
     +RR,GofR,CONU1,CONU2,RKBI,'rdf:',NM(IS),'-',NM(JS),IRDF(K,I)
            end if
	    if(I.lt.IFST)IFST=I
	   end if
	 end do
	 if(NPRI.ge.1.and..not.(LCOMRDF1.and.LCOMRDF2))then
	   write(16,*)
	   write(16,*)' Intramolecular RDF   '
	   write(16,*)' ----------------------------------------- '
           write(16,'(a)')
     + '  R          RDF       dist.distribution'
	write(16,*)'-------------------------------------------------'
 	   FACTI=1.d0/(NSPEC(ITS(IS))*NPRI*IAN)
	   do I=1,NAI
	     RR=RMI+(I-0.5)*DRI/NAI
	     DRR=DRI/NAI
   	     SH12=4*PI*DRR*(RR**2+DRR**2/12.d0)
	     DIST=IRDFI(K,I)*FACTI/DRR
             GofR=IRDFI(K,I)*FACT/SH12
	     if(DIST.gt.1.d-8)then
	      write(16,'(f8.3,2f13.5,11x,a8,2x,a4,a1,a4)')
     +     RR,GofR,DIST,'rdf_int:',NM(IS),'-',NM(JS)
	     end if
	   end do
	 end if
      end do
      write(16,*)' SUMMARY:'
      write(16,'(a1,7x,50(1x,a4,a1,a4))')
     +  '#',(NM(IRD(K,1)),'-',NM(JRD(K,1)),K=1,NRDF)
      do I=IFST,NA
	RR=(I-0.5)*RDFCUT/NA
	DRR=RDFCUT/NA
  	SH12=4*PI*DRR*(RR**2+DRR**2/12.d0)
        do K=1,NRDF
	 NPAIR=0
	 do KK=1,IREP(K)
	   IS=IRD(K,KK)
	   JS=JRD(K,KK)
	   if(IS.eq.JS)then
	     NPAIR=NPAIR+NSPEC(ITS(IS))*(NSPEC(ITS(JS))-1)/2
	   else
	     NPAIR=NPAIR+NSPEC(ITS(IS))*NSPEC(ITS(JS))
	   end if
	 end do
	 FACT=VOLS/(NPAIR*1.d0*IAN)
	 TMP(K)=IRDF(K,I)*FACT/SH12
	end do
	write(16,'(f8.3,50f10.5)')RR,(TMP(K),K=1,NRDF)
      end do
      close(16)
      stop
      end
