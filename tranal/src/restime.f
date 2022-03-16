      program RESTIME
*  MDynaMix trajectory analysis suite
*  From version 5.2
*  Calculation of residence time
      include "tranal.h"
C     array boundaries for:
C     NRT1M: number of sites "from"
C     NRT2M: numberof sites "to"
C     NITTM: number of time frames to follow residence time      
      parameter (NRT1M=4096,NRT2M=8192)
      real*8 TIIN(NRT1M,NRT2M),TIINI(NRT1M,NRT2M),TIINO(NRT1M,NRT2M)
      real*8, allocatable :: DISTT(:),DIST2(:)
      integer IRS(NRT1M,NRT2M),IR1(NRT1M),IR2(NRT2M),IRS1(NRT1M),
     & IRS2(NRT2M)
      character*64 FILOUT
      namelist /RESTIM/NR1,NR2,IR1,IR2,NITT,RRN,RRX,RRD,TTER,TIMM,FILOUT
      data IINIT/0/
*  Trajectory data
      call SETTRAJ
*  Default parameters
      TTER = 2.       ! 2 ps
      RRD = 0.25
      RRN = 0.
*  Input of specific parameters
      read(*,RESTIM)
      allocate(DISTT(NITT),DIST2(NITT))
      if(RRN.gt.RRX) stop 'RRN>RRX'
      RRNN=RRN-RRD
      RRNX=RRN+RRD
      RRXN=RRX-RRD
      RRXX=RRX+RRD
*  Compose list of atoms for calculation of residence time
	NRT1=0
	NRT2=0
	IHP=0
	do I=1,NR1
	  IS=IR1(I)
	  ITYP1=ITS(IS)
	  NS1=NSPEC(ITYP1)
	  do J=1,NS1
	    NRT1=NRT1+1
	    if(NRT1.gt.NRT1M)stop ' Increase NRT1M  in restime.f'
	    ISP=ISADDR(ITYP1)+(J-1)*NSITS(ITYP1)+IS-ISADR(ITYP1)
	    IRS1(NRT1)=ISP
	  end do
      end do
      do I=1,NR2
	IS=IR2(I)
	ITYP2=ITS(IS)
	NS2=NSPEC(ITYP2)
	do J=1,NS2
          NRT2=NRT2+1
	  if(NRT2.gt.NRT2M)stop ' increase NRT2M in restime.f'
	  ISP=ISADDR(ITYP2)+(J-1)*NSITS(ITYP2)+IS-ISADR(ITYP2)
	  IRS2(NRT2)=ISP
	end do
      end do
      if(IPRINT.ge.8)write(*,*)NRT1,(IRS1(I),I=1,NRT1)
      if(IPRINT.ge.8)write(*,*)NRT2,(IRS2(I),I=1,NRT2)
*  Starting parameters
	IEND = 0
      IBREAK = 0
      IMD=0
      ITIR=0
      ITIR1=0
      ITIR2=0
      TIMRS=0.
      TIMRS1=0.
      TIMRS2=0.
      do I=1,NITTM
	DISTT(I)=0.
	DIST2(I)=0.
      end do
*  Begin to read trajectory
      do while(IEND.eq.0)
        call READCONF(IEND)
	if(IEND.ne.0)IMD=1
	if(IINIT.eq.0)then
	  IINIT=1
	  do I=1,NRT1M
	    do J=1,NRT2M
	      IRS(I,J)=0
	    end do
	  end do 
	end if 
	if(IBREAK.ne.0)then
	  do I=1,NRT1
	    do J=1,NRT2
	      IRS(I,J)=0
	    end do
	  end do 
	  write(*,*)' residence time tables set to zero'
	end if     
*   over first set of sites
	do I=1,NRT1
	  ISP=IRS1(I)
	  ITYP=ITYPE(ISP)
	  IMOL=NNUM(ISP)
	  if(LIST(ITYP).eq.0)go to 99
*   over second set of sites
	  do J=1,NRT2
	    JSP=IRS2(J) 
	    JMOL=NNUM(JSP)
	    if(IMOL.ne.JMOL)then
	    JTYP=ITYPE(JSP)
	    if(LIST(JTYP).eq.0)go to 99
	    if(ISP.eq.JSP)go to 90
	    DX=SX(ISP)-SX(JSP)
	    DY=SY(ISP)-SY(JSP)
	    DZ=SZ(ISP)-SZ(JSP)
	    call PBC(DX,DY,DZ)
	    RR=sqrt(DX**2+DY**2+DZ**2)
*  executed at the last call
	    if(IMD.eq.1)then
*  register all pairs which are still in contact
	      if(mod(IRS(I,J),3).ne.0)then
	        TIMR=FULTIM-TIINI(I,J)
	        TIMRS2=TIMRS2+TIMR     ! register residence time 
	        NR=TIMR*NITT/TIMM+1
		write(*,*)' add ',TIMR,' ps  for ',I,J,' R=',RR
	        if(NR.gt.0.and.NR.le.NITT)DIST2(NR)=DIST2(NR)+1 
	        ITIR2=ITIR2+1
	        if(TIMR.gt.TTER)then
	          TIMRS1=TIMRS1+TIMR
	 	  ITIR1=ITIR1+1
		end if
	        TIMRS=TIMRS+TIMR     
	        if(NR.gt.0.and.NR.le.NITT)DISTT(NR)=DISTT(NR)+1 
	        ITIR=ITIR+1
	      end if
	      go to 90
	    end if    ! IMD.eq.1
	    if(IPRINT.ge.8)write(*,*)I,J,ISP,JSP,RR
*  standard criteria (exit during short time not accounted)
 15	    if(mod(IRS(I,J),3).eq.0)then 
*  bulk (0) -> enter hydrated shell (2)
	      if(RR.gt.RRN.and.RR.lt.RRX)then
		IRS(I,J)=IRS(I,J)+2
	 	TIINI(I,J)=FULTIM          ! register entrance time
	        if(IPRINT.ge.7)write(*,*)' new st',I,J
	      end if   
*  Temporarly out (1)
	    else if(mod(IRS(I,J),3).eq.1)then 
*  Out more than cut-off time (TTER) -> bulk (0) 
	      if(FULTIM-TIINO(I,J).gt.TTER)then
	        IRS(I,J)=IRS(I,J)-1
	        if(IPRINT.ge.7)write(*,*)' out st',I,J
	        TIMR=(TIINO(I,J)-TIINI(I,J))
	        TIMRS2=TIMRS2+TIMR     ! register residence time 
	        NR=TIMR*NITT/TIMM+1
	        if(NR.gt.0.and.NR.le.NITT)DIST2(NR)=DIST2(NR)+1 
	        ITIR2=ITIR2+1
	        if(TIMR.gt.TTER)then
	          TIMRS1=TIMRS1+TIMR
	 	  ITIR1=ITIR1+1
		end if
	 	go to 15
	      end if
	      if(RR.gt.RRN.and.RR.lt.RRX)IRS(I,J)=IRS(I,J)+1
*  hydrated (2)
	    else
	      if(RR.le.RRN.or.RR.ge.RRX)then
*  mark as temporarly out (1); register exit time 
	        IRS(I,J)=IRS(I,J)-1
	        TIINO(I,J)=FULTIM
              end if
	    end if	    
*  alt. criteria (just register in - out )	     
	    if(IRS(I,J).le.2)then
* simple in
	      if(RR.gt.RRNX.and.RR.lt.RRXN)then
	        IRS(I,J)=IRS(I,J)+3
	        TIIN(I,J)=FULTIM 
	        if(IPRINT.ge.7)write(*,*)' new alt',I,J
	      end if
	    else
* simple out    
	      if(RR.lt.RRNN.or.RR.gt.RRXX)then 
	        TIMR=(FULTIM-TIIN(I,J))
		IRS(I,J)=IRS(I,J)-3
	        TIMRS=TIMRS+TIMR     
	        NR=TIMR*NITT/TIMM+1
	        if(NR.gt.0.and.NR.le.NITT)DISTT(NR)=DISTT(NR)+1 
	        ITIR=ITIR+1
	        if(IPRINT.ge.7)write(*,*)' out alt',I,J
	      end if
	    end if
*  Count hydrated pairs
	    if(mod(IRS(I,J),3).eq.2)IHP=IHP+1
 90	    continue
	  end if ! (IMOL.ne.JMOL
	  end do
	end do
 99	if(IPRINT.ge.6)write(*,*)fultim,ITIR,TIMRS/ITIR
      end do
C  begin output
 100  continue
      if(IR1(1).eq.IR2(1))IHP=IHP/2
      open(unit=15,file=filout,status='unknown')
      write(15,'(a)')'#  Distribution over residence times'
      write(15,'(a,I5,a,I5)')'#  sites ',IR1(1),' and ',IR2(1)
      write(15,'(a,f9.3,a,f9.3,a,f7.3)')
     & '#  between ',RRN,' and ',RRX,' +/- ',RRD 
      write(15,'(a,f8.2,a)')'#  permission time T* = ',TTER,' ps '
      write(15,'(a,I7,a,i7,a)')'#  total ',ITIR2,' and ',ITIR,' events'
      write(15,'(a,I7)')'#  cases of longer T* time ',ITIR1
      write(15,'(a,f10.3)')'#  average number of bound pairs ',
     +           IHP*1./IAN
      write(15,'(a,f10.3)')'#  bound pairs per site 1:       ',
     +           IHP*1./(NS1*IAN)
      write(15,*)
      write(15,'(a,f14.3,a)')
     + '#  Average residence time (with permission) ',
     +          TIMRS2/ITIR2,' ps' 
      write(15,'(a,f14.3,a)')
     +  '#  Average residence time (simple average)  ',TIMRS/ITIR,' ps' 
      write(15,'(a,f14.3,a)')
     +  '#  Exponential decay time                   ',
     +  TIMRS1/ITIR1-TTER,' ps'
	  write(15,'(a)')
     +'                With permission                  Simple av.'
	  write(15,'(a)')
     +'#  time(ps)     cases-1      P1      P(remain)     cases-2   P2'
	  DTIRS=ITIR2
	  do I=1,NITT  
	    PPP1=DIST2(I)/ITIR2
	    DTIRS=DTIRS-DIST2(I)
	    PFUN=DTIRS/ITIR2          
	    TTT=I*TIMM/NITT                 
	    PPP=DISTT(I)/ITIR          
	    write(15,'(f12.3,2x,f7.0,2f13.8,2x,f7.0,f13.8)')
     +    TTT,DIST2(I),PPP1,PFUN,DISTT(I),PPP
	  end do
	stop
	end
