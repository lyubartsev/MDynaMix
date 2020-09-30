	program MAKEMOL              
C
C   makemol v. 5.2.7
C
C   This utility allows to construct .mmol files for the MDynaMix program
C   from combination of two files: simple molecular file (.smol)
C   describing only molecular coordinates, atom names and list of bonds,
C   and second file describing the force field. For additional information, 
C   see Makemol.doc file.
C
C   These parameters may be changed:
	parameter (NBB=20000,NAA=30000,NTT=50000,NNBM=200,NAT=100000,
     + NP14M=50000)
C
	real*8 RB(NBB),FB(NBB),RA(NAA),FA(NAA),SA(NAA),FSA(NAA),
     +RT(NTT),FT(NTT),SGM(NNBM),EPS(NNBM),SGM14(NNBM),EPS14(NNBM),
     +DB(NBB),R0(NBB),SGP14(NP14M),EPP14(NP14M),MASS(NAT)
	real*8 X(NAT),Y(NAT),Z(NAT),SM(NAT),Q(NAT),EP(NAT),SI(NAT)
	character*4 NB1(NBB),NB2(NBB),NA1(NAA),NA2(NAA),NA3(NAA),
     +NT1(NTT),NT2(NTT),NT3(NTT),NT4(NTT),NPI(NP14M),NPJ(NP14M)
	character*4 NM(NAT),NFF(NAT),NN(NAT),NFM(NAT),
     +NDUB,CDUB*8,CDUB2*8
	integer NMUL(NTT),IER(NAT),IB(NBB),JB(NBB),ITMP(20),
     +IBUK(NBB),IAUK(NAA),ITMP2(20),ITUK(NTT),ITUKA(NTT),
     +IA1(NAA),IA2(NAA),IA3(NAA),IT1(NTT),IT2(NTT),IT3(NTT),IT4(NTT),
     +IUK14(NAT),IUKUB(NAA),IREF(NAT),IU14(NP14M),JU14(NP14M),
     +MU14(NP14M)
	character*128 FIN,FDB,FOUT,STR,AUX,TAKESTR
*
*    1. Input file name information
*
	numargs=iargc()
	if(numargs.ne.2)then
	   write(*,*)'usage:'
	   write(*,*)'makemol <mol file> <force field file>'
	   stop
	end if
	call getarg(1,FIN)
	call getarg(2,FDB)
	LI=LENG(FIN,128)
	LF=LENG(FDB,128)
	if(LI.le.4)then
	   write(*,*)'Molecular structure file type not recognized'
	   write(*,*)'File name ',FIN(1:LI)
	   stop
	end if
	if(FIN(LI-4:LI).eq.'.smol')then
	   IINP=0
	   write(*,*)'Assuming .smol molecular structure file format'
	   FOUT(1:LI-5)=FIN(1:LI-5)
	   FOUT(LI-4:LI)='.mmol'
	else if(FIN(LI-3:LI).eq.'.crd')then
	   IINP=1
	   write(*,*)'Assuming .crd molecular structure file format'
	   FOUT(1:LI-4)=FIN(1:LI-4)
	   FOUT(LI-3:LI+1)='.mmol'
	else if(FIN(LI-3:LI).eq.'.alc')then
	   IINP=2
	   write(*,*)'Assuming .alc molecular structure file format'
	   FOUT(1:LI-4)=FIN(1:LI-4)
	   FOUT(LI-3:LI+1)='.mmol'
	else if(FIN(LI-4:LI).eq.'.mol2')then
	   IINP=3
	   write(*,*)'Assuming .mol2 molecular structure file format'
	   FOUT(1:LI-5)=FIN(1:LI-5)
	   FOUT(LI-4:LI)='.mmol'
	else
	   write(*,*)'Molecular structure file type not recognized'
	   write(*,*)'File name ',FIN(1:LI)
	   write(*,*)'Acceptable formats: .mol2 .smol .crd .alc'
	   stop
	end if
	open(unit=20,file='makemol.log',status='unknown')
	open(unit=11,file=FIN,status='old',err=11)
	go to 12
 11	write(*,*)' Molecular structure file ',FIN,' not found'
 	write(20,*)' Molecular structure file ',FIN,' not found'
	stop
 12	open(unit=12,file=FDB,status='old',err=14)
 	go to 15
 14	write(*,*)' Force field file ',FDB,' not found'
 	write(20,*)' Force field file ',FDB,' not found'
	stop
 15	continue
	write(20,*)' Open Force Field file ',FDB(1:LF)
*
*    2. Read force field database
*    ----------------------------
	NBD=0 
	NAN=0
	NTOR=0
	NBN=0
	NP14=0
	fudgeLJ=1.
	fudgeQQ=1.
	NSTS=0
	ISKIP=0
	IUB=1
	IAM=0
*
*    2.1 Bonds
*
 20	STR=TAKESTR(12,IE)
 	if(STR(1:9).eq.'begin_def')then
	   IA=1
 23	   STR=TAKESTR(12,IE)
	   if(STR(1:5).eq.'No_UB')IUB=0
	   if(STR(1:8).eq.'SKIP_ABS')ISKIP=1
           if(STR(1:7).eq.'FACTORS')then
	      read(STR(9:80),*)fudgeLJ,fudgeQQ
	      write(*,'(a,2f8.4)')
     &' 1-4 general scaling factors: ',fudgeLJ,fudgeQQ
	      write(20,'(a,2f8.4)')
     &' 1-4 general scaling factors: ',fudgeLJ,fudgeQQ
	   end if 
	   if(STR(1:7).eq.'end_def')go to 21
	   read(STR,*,err=23,end=23)NFM(IA),MASS(IA)
	   IAM=IA
	   IA=IA+1
           go to 23
	end if
 21	if(STR(1:4).eq.'BOND')then
 22	  I=NBD+1
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=26,end=26)NB1(I),NB2(I),RB(I),FB(I),DB(I),R0(I)
	  go to 27
 26	  read(STR,*,err=28,end=28)NB1(I),NB2(I),RB(I),FB(I)
	  DB(I)=0.
	  R0(I)=0.
 27	  NBD=NBD+1
	  if(NBD.gt.NBB)stop 'Increase NBB'
	  go to 22
 28	  write(*,*)' read ',NBD,' ff bonds'
 	  write(20,*)' read ',NBD,' ff bonds'
*
*    2.2 Angles
*
	else if(STR(1:4).eq.'ANGL')then
 31	   I = NAN+1
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=36,end=36)
     +     NA1(I),NA2(I),NA3(I),RA(I),FA(I),SA(I),FSA(I)
	  go to 37 
 36	  read(STR,*,err=38,end=38)NA1(I),NA2(I),NA3(I),RA(I),FA(I)     
	  SA(I)=-1.
	  FSA(I)=0.
 37	  NAN = NAN+1
	  if(NAN.gt.NAA)stop 'Increase NAA'
          go to 31
 38	  write(*,*)' read ',NAN,' ff angles'
 	  write(20,*)' read ',NAN,' ff angles'
*	  
*    2.3 Torsions
*
	else if(STR(1:4).eq.'TORS'.or.STR(1:4).eq.'DIHE')then
 41	   I=NTOR+1
	  STR=TAKESTR(12,IE)
	  read(STR,*,end=48,err=48)
     +    NT1(I),NT2(I),NT3(I),NT4(I),RT(I),FT(I),NMUL(I)
	  NTOR=NTOR+1
	  if(NTOR.gt.NTT)stop 'Increase NTT'
	  go to 41
 48	  write(*,*)' read ',NTOR,' ff torsions'
 	  write(20,*)' read ',NTOR,' ff torsions'
*
*    2.4 Non-bonded
*
	else if(STR(1:4).eq.'NONB')then
 51	  I=NBN+1
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=56,end=56)
     +     NN(I),SGM(I),EPS(I),SGM14(I),EPS14(I)
	  go to 57 
 56	  read(STR,*,err=58,end=58)NN(I),SGM(I),EPS(I)
	  if(fudgeLJ.eq.1.)then
 	    SGM14(I)=-1.
	    EPS14(I)=-1.
	  else
	    SGM14(I)=SGM(I)
	    EPS14(I)=fudgeLJ*EPS(I)
	  end if
 57	  NBN=NBN+1
	  if(NBN.gt.NNBM)stop 'Increase NNBM'
	  go to 51
 58	  write(*,*)' read ',NBN,' ff non-bonded'
 	  write(20,*)' read ',NBN,' ff non-bonded'
*
*    2.5  1-4 pairs
*
	else if(STR(1:5).eq.'PAR14')then
 73	   I=NP14+1
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=76,end=76)
     +     NPI(I),NPJ(I),SGP14(I),EPP14(I)
	  NP14=NP14+1
	  if(NP14.gt.NP14M)stop 'Increase NP14M'
	  go to 73
 76	  write(*,*)' read ',NP14,' 1-4 pairs'
	  write(20,*)' read ',NP14,' 1-4 pairs'	
	else if(STR(1:5).eq.'LJ_PAIRS')then
 83	  STR=TAKESTR(12,IE)
 	  read(STR,*,err=86,end=86)CDUB,CDUB2,SIG,EPS
	  go to 83
 86	  write(*,*)' Read non-nobded LJ pair section'
	  write(20,*)' Read non-nobded LJ pair section'
*        
*  2.6 Finish reading database
*
	else if(STR(1:3).eq.'END')then
	  write(*,*)' FF parameters from database ',FDB,' loaded'
	  go to 60
	else
	  go to 20
	end if
	go to 21
 60	close(12)
	write(20,*)' '
*
	open(unit=13,file=FOUT,status='unknown') 
	write(*,*)' Writing molecular topology-FF file ',FOUT
	write(20,*)' Writing molecular topology-FF file ',FOUT
*
*     3. Read simple molecular file
*
	if(IINP.eq.3)then
*     3.1   mol2.  format
 61	  read(11,'(a100)',end=85,err=99)AUX
          if(AUX(1:1).ne.'#'.and.AUX(1:1).ne.'*')go to 62
          LN=LENG(AUX,128)
	  AUX(1:1)='#'
          write(13,'(a)')AUX(1:LN)
          go to 61
 62	  if(AUX(1:9).eq.'@<TRIPOS>')then
	    if(AUX(10:13).eq.'MOLE')then
	      read(11,*)AUX
              LN=LENG(AUX,128)
              write(13,'(a,a)')'# ',AUX(1:LN)
	      read(11,*,err=199,end=199)NSTS,NB
	      if(NB.gt.NBB)stop 'increase NBB'
	      write(*,*)'  MOLECULES Section is read '
	      write(*,'(2(i6,a))')NSTS, ' bonds, ',NB, ' bonds' 
	      write(20,*)'  MOLECULES Section is read '
	      write(20,'(2(i6,a))')NSTS, ' bonds, ',NB, ' bonds' 
	    end if
	    if(AUX(10:13).eq.'ATOM')then
	      do I=1,NSTS
	read(11,*,err=63,end=63)J,NM(I),X(I),Y(I),Z(I),NFF(I),II,NDUB,Q(I)
*		write(*,*)J,NM(I),X(I),Y(I),Z(I),NFF(I),II,NDUB,Q(I)
		IREF(J)=I
              end do
	      go to 64
 63	      write(*,*)' Error in the list of atoms, atom nr ',I
	      stop
 64	      continue
	      write(*,*)'  ATOMS Section is read '
	      write(20,*)'  ATOMS Section is read '
	    end if
	    if(AUX(10:13).eq.'BOND')then
	      do I=1,NB
	        read(11,*,end=299,err=299)J,IBB,JBB
	        IB(I)=IREF(IBB)
	        JB(I)=IREF(JBB)
	        if(IB(I).le.0.or.IB(I).gt.NSTS.or.JB(I).le.0.or.
     +    JB(I).gt.NSTS.or.IB(I).eq.JB(I))then
	          write(*,*)' Bad bond!!! ',I,IB(I),JB(I)
	          write(20,*)' Bad bond!!! ',I,IB(I),JB(I)
	        end if
	      end do
	      write(*,*)'  BONDS Section is read '
	      write(20,*)'  BONDS Section is read '
	    end if
	  end if
	  go to 61
	else
*     3.2.This keep initial commentary lines
 70       read(11,'(a100)',end=99,err=99)AUX
          if(AUX(1:1).ne.'#'.and.AUX(1:1).ne.'*')go to 71
          LN=LENG(AUX,128)
	  AUX(1:1)='#'
          write(13,'(a)')AUX(1:LN)
          go to 70
*      3.2.2 Number of atoms (first uncommented line)
 71       read(AUX,*)NSTS      
	  write(*,*)' num of sites ',NSTS
	  write(20,*)' num of sites ',NSTS
*      3.2.3 Atom names and coordinates
	  do I=1,NSTS
	    IER(I)=0
            AUX=TAKESTR(11,IE)
	    LN=LENG(AUX,100)
	    if(IINP.eq.1)then
	      read(AUX,*,end=199,err=199)
     +        J,ND1,NDUB,NM(I),X(I),Y(I),Z(I),NFF(I),ND3,Q(I)
	      IREF(J)=I
	    else if(IINP.eq.2)then
	      read(AUX,*,end=199,err=199)
     +        J,NFF(I),X(I),Y(I),Z(I),Q(I),NM(I)
	      IREF(J)=I
	    else
	      read(AUX,*,end=199,err=199)NM(I),X(I),Y(I),Z(I),Q(I),NFF(I)
	      IREF(I)=I
	    end if
	  end do
*      3.2.4 Bond list
	  AUX=TAKESTR(11,IE)
	  read(AUX,*,end=199,err=199)NB
	  if(NB.gt.NBB)stop 'increase NBB'
	  IOK=0
	  do I=1,NB
	    AUX=TAKESTR(11,IE)
	    IOK=IOK+1
	    if(IINP.eq.2)then
	      read(AUX,*,end=199,err=199)J,IBB,JBB
	      IB(I)=IREF(IBB)
	      JB(I)=IREF(JBB)
	    else
	      read(AUX,*,end=199,err=199)IBB,JBB
	      IB(I)=IREF(IBB)
	      JB(I)=IREF(JBB)
	    end if
	    if(IB(I).le.0.or.IB(I).gt.NSTS.or.JB(I).le.0.or.
     +    JB(I).gt.NSTS.or.IB(I).eq.JB(I))then
	     write(*,*)' Bad bond!!! ',I,IB(I),JB(I)
	     write(20,*)' Bad bond!!! ',I,IB(I),JB(I)
	     IOK=IOK-1
	    end if
	  end do
	  NB=IOK
	end if
 85	NUB=0
*  
*    4. Creating new .mmol file
*    --------------------------
*    4.1 Finding information on masses, sigma, epsilon
*    -------------------------------------------------
*    4.1.1. Masses
	do I=1,NSTS
C  First looking fo the definitions
	  do J=1,IAM
	    if(NFM(J).eq.NM(I))then
	      SM(I)=MASS(J)
	      go to 79
	    end if
	  end do
C  Now guessing from the names
	  if(NM(I)(1:1).eq.'H')then
	     if(NM(I)(2:2).eq.'e'.or.NM(I)(2:2).eq.'E')then
C  Helium		
	       SM(I)=4.
	     else
C  Hydrogen
	       SM(I)=1.008
	     end if
	  else if(NM(I)(1:2).eq.'Al'.or.NM(I)(1:2).eq.'AL')then
C  Alluminium
	     SM(I)=26.98
	  else if(NM(I)(1:1).eq.'O')then
C  Oxygen
	     SM(I)=15.9994
	  else if(NM(I)(1:2).eq.'Br'.or.NM(I)(1:2).eq.'BR'.or.
     +            NM(I)(1:2).eq.'br')then
C   Brom
	     SM(I)=79.9
	  else if(NM(I)(1:1).eq.'C')then
C  Carbon
	     SM(I)=12.011
	  else if(NM(I)(1:1).eq.'F')then
	     if(NM(I)(1:2).eq.'FE'.or.NM(I)(1:2).eq.'Fe')then
C  Iron
	       SM(I)=55.847
	     else
C  Floor
	       SM(I)=19.
	     end if
	  else if(NM(I)(1:1).eq.'N')then
C  Nitrogen
	     SM(I)=14.007
	  else if(NM(I)(1:1).eq.'P')then
C  Phosphorus
	     SM(I)=30.974
	  else if(NM(I)(1:1).eq.'S')then
	     if(NM(I)(2:2).eq.'i'.or.NM(I)(2:2).eq.'I')then
C  Silicon		
	       SM(I)=28.09
	     else
C  Sulfur
	       SM(I)=32.06
	     end if
	  else if(NM(I)(1:2).eq.'TI'.or.NM(I)(1:2).eq.'Ti')then
C  Titanium
	     SM(I)=47.87
	  else if(NM(I)(1:2).eq.'RU'.or.NM(I)(1:2).eq.'Ru')then
C  Ruthenuim
	     SM(I)=101.07
	  else 
C  Unknown atom
	     write(*,*)' Unknown atom ',NM(I),' num ',I
	     write(20,*)' Unknown atom ',NM(I),' num ',I
	     SM(I)=10.
*	     IER(I)=1
	  end if
 79	  continue
	end do
*  4.1.2 Sigma and Epsilon
	write(13,'(a)')'#  Number of atoms '
	write(13,*)NSTS
	NSP14=0
	do I=1,NSTS
	  do J=1,NBN
	     if(NFF(I).eq.NN(J))go to 81
	  end do
	  write(*,'(a,i6,2x,a)')
     +' Nonbonded parameters not found for atom ',I,NFF(I)
	  write(20,'(a,i6,2x,a)')
     +' Nonbonded parameters not found for atom ',I,NFF(I)
	  IER(I)=1
	  EP(I)=0.
	  SI(I)=0.
	  go to 82
 81	  EP(I)=EPS(J)
	  SI(I)=SGM(J)
	  if(EPS14(J).gt.0..and.EPS14(J).gt.0.)then
	    NSP14=NSP14+1
	    IUK14(I)=J
	  else
	    IUK14(I)=0
	  end if
C  writing to .mmol file
 82	  if(IER(I).eq.0)then
	    write(13,
     +'(a4,3(1x,f8.3),4(1x,f8.4),1x,a4,1x,i5)')
     +      NM(I),X(I),Y(I),Z(I),SM(I),Q(I),SI(I),EP(I),NFF(I),I
	  else
	    write(13,
     +'(a1,a4,3(1x,f8.3),4(1x,f8.4),1x,a4,1x,i5)')
     +      '!',NM(I),X(I),Y(I),Z(I),SM(I),Q(I),SI(I),EP(I),NFF(I),I
	  end if
	end do
	write(*,*)' found ',NSP14,' atoms with special 1-4 parameters' 
*  4.1.3  Commentary
	write(13,*)'3'
	write(13,*)'File created by utility MAKEMOL, package MDynaMix'
	write(13,'(a,a)')'Coordinates taken from file ',FIN(1:LI)
	write(13,'(a,a)')'Force field from file ',FDB(1:LF)
*
*   4.2. Bonds
*   ----------
	do I=1,NB
	  do J=1,NBD
	    IBI=IB(I)
	    JBJ=JB(I)
	    if(NFF(IBI).eq.NB1(J).and.NFF(JBJ).eq.NB2(J).or.
     +       NFF(JBJ).eq.NB1(J).and.NFF(IBI).eq.NB2(J))go to 95
	  end do
	  write(*,*)' parameters for bond ',NM(IBI),'-',NM(JBJ),
     +' or ',NFF(IBI),'-',NFF(JBJ),' not found'
	  write(20,*)' parameters for bond ',NM(IBI),'-',NM(JBJ),
     +' or ',NFF(IBI),'-',NFF(JBJ),' not found'
	  IBUK(I)=0
	  go to 96
 95	  IBUK(I)=J
 96	  continue
	end do
*
*    4.3 Angles
*    ----------
*    4.3.1 Compose list of anglefrom the list of bonds
C    we go through the list of atoms and consider all the bonds
	NA=0
	NAP=0
	do I=1,NSTS
	   NBONDS=0
	   do J=1,NB
	     if(IB(J).eq.I)then
		NBONDS=NBONDS+1
		ITMP(NBONDS)=JB(J)
             end if
	     if(JB(J).eq.I)then
		NBONDS=NBONDS+1
		ITMP(NBONDS)=IB(J)
             end if
	   end do
C  this atom has NBONDS bound atoms
	   if(NBONDS.ge.2)then
	      do J1=1,NBONDS-1
		 do J2=J1+1,NBONDS
		    NA=NA+1
		    IA1(NA)=ITMP(J1)
		    IA2(NA)=I
		    IA3(NA)=ITMP(J2)
		 end do
	      end do
	   end if
	end do
	write(*,*)' Total number of angles ',NA
	write(20,*)' Total number of angles ',NA
	if(NA.gt.NAA)stop 'increase NAA'
*
*   4.3.2 Finding parameters for angles
*
	do I=1,NA
	  IER(I)=0
	  IA=IA1(I)
	  JA=IA2(I)
	  KA=IA3(I)
	  do J=1,NAN
	    if((NFF(IA).eq.NA1(J).and.NFF(JA).eq.NA2(J).and.
     +        NFF(KA).eq.NA3(J)).or.(NFF(KA).eq.NA1(J).and.
     +        NFF(JA).eq.NA2(J).and.NFF(IA).eq.NA3(J)))go to 105
	  end do
	  write(*,*)' parameters for angle ',NM(IA),'-',NM(JA),'-',
     +NM(KA),' or ',NFF(IA),'-',NFF(JA),'-',NFF(KA),' not found'
	  write(20,*)' parameters for angle ',NM(IA),'-',NM(JA),'-',
     +NM(KA),' or ',NFF(IA),'-',NFF(JA),'-',NFF(KA),' not found'
	  IAUK(I)=0
	  IUKUB(I)=0
	  go to 106
 105	  IAUK(I)=J
	  NAP=NAP+1
	  if(SA(J).gt.0.and.IUB.gt.0)then
C
C  Urey-Bradley term (harmonic bond for 1-3 interactions)
C
	     IUKUB(I)=J
	     NUB=NUB+1
	  else
	     IUKUB(I)=0
	  end if
 106	  continue
        end do
	write(20,*)' parameters for ',NAP, 'angles found'
	if(ISKIP.eq.1)write(*,*)' parameters for ',NAP, 'angles found'
	if(IUB.gt.0)write(20,*)NUB,' Urey-Bradley terms added'
*
*   4.4 Add bonds to .mmol file
*   --------------------------
      write(13,'(a)')'#  Number of bonds '
      write(13,*)NB+NUB
      do I=1,NB
	J=IBUK(I)
	IBI=IB(I)
	JBJ=JB(I)
	if(J.gt.0)then
	  RRB=RB(J)
	  FFB=FB(J)
	  if(DB(J).le.1.d-6)then
	    JD=0 
  	    write(13,'(3I6,f9.5,f10.2,21x,a1,2x,a4,a1,a4)') 
     +JD,IB(I),JB(I),RRB,FFB,'#',NFF(IBI),'-',NFF(JBJ)
	  else
	     JD=1
             DBB=DB(J)
	     RO=R0(J)
  	write(13,'(3I6,f9.5,f10.2,1x,2f9.4,2x,a1,2x,a4,a1,a4)')
     + JD,IB(I),JB(I),RRB,FFB,DBB,RO,'#',NFF(IBI),'-',NFF(JBJ) 
	  end if
	else
	  JD=0
	  RRB=0.
	  FFB=0.
          write(13,'(a1,3I6,f9.5,f10.2,21x,a1,2x,a4,a1,a4)') 
     +'!',JD,IB(I),JB(I),RRB,FFB,'#',NFF(IBI),'-',NFF(JBJ)
	end if
      end do 
      IAD=0
*  Adding Urey-Bradley terms
      if(IUB.gt.0)then
	write(13,'(a)')'#    Urey-Bradley terms (for 1-3 neigbours)'
        do I=1,NA
	  J=IUKUB(I)
	  if(J.gt.0)then
	    IAD=IAD+1
	    JD=2
	    RRB=SA(J)
	    FFB=FSA(J)
  	    write(13,'(3I6,f9.5,f10.2,21x,a1,2x,a4,a1,a4)') 
     +JD,IA1(I),IA3(I),RRB,FFB,'#',NFF(IA1(I)),'-',NFF(IA3(I))
	  end if
        end do
      end if
      if(IAD.ne.NUB)then
	 write(*,*)' Error: Number of UB terms wrong:',IAD,NUB
	 write(20,*)' Error: Number of UB terms wrong:',IAD,NUB
      end if
      write(*,*)' bond list completed'
      write(20,*)' bond list completed'
*
*   4.5 Add angles to .mmol file
*   ---------------------------
      write(13,'(a)')'#  Number of angles '
      if(ISKIP.eq.1)then
	write(13,*)NAP
      else
        write(13,*)NA
      end if
      do I=1,NA
	J=IAUK(I)
	IAI=IA1(I)
	JAJ=IA2(I)
	KAK=IA3(I)
	if(J.gt.0)then
	  RRA=RA(J)
	  FFA=FA(J)
  	  write(13,'(3I6,f10.4,f10.2,5x,a1,2x,3(a4,a1))') 
     +IAI,JAJ,KAK,RRA,FFA,'#',NFF(IAI),'-',NFF(JAJ),'-',NFF(KAK)
	else
	  RRA=0.
	  FFA=0.
  	  if(ISKIP.ne.1)
     +write(13,'(a1,3I6,f10.4,f10.2,5x,a1,2x,3(a4,a1))') 
     +'!',IAI,JAJ,KAK,RRA,FFA,'#',NFF(IAI),'-',NFF(JAJ),'-',NFF(KAK)
	end if
      end do                    
      write(*,*)' angles list completed'
      write(20,*)' angles list completed'

*
*   4.6 Torsions
*   -----------
*   4.6.1 Prepare list of torsions
	NT=0
	NTP=0
	do I=1,NB
	   NBD1=0
	   NBD2=0
	   I2=IB(I)
	   I3=JB(I)
	   do J=1,NB
	     if(IB(J).eq.I2)then
	       if(JB(J).ne.I3)then
		 NBD1=NBD1+1
		 ITMP(NBD1)=JB(J)
	       end if
	     end if
	     if(IB(J).eq.I3)then
	       if(JB(J).ne.I2)then
		 NBD2=NBD2+1
		 ITMP2(NBD2)=JB(J)
	       end if
	     end if
	     if(JB(J).eq.I2)then
	       if(IB(J).ne.I3)then
		 NBD1=NBD1+1
		 ITMP(NBD1)=IB(J)
	       end if
	     end if
	     if(JB(J).eq.I3)then
	       if(IB(J).ne.I2)then
		 NBD2=NBD2+1
		 ITMP2(NBD2)=IB(J)
	       end if
	     end if
	   end do
	   if(NBD1*NBD2.gt.0)then
	      do J=1,NBD1
		do K=1,NBD2
		   NT=NT+1
		   IT1(NT)=ITMP(J)
		   IT2(NT)=I2
		   IT3(NT)=I3
		   IT4(NT)=ITMP2(K)
		end do
              end do
	   end if
	end do
	write(*,*)' Total number of torsions ',NT
	write(20,*)' Total number of torsions ',NT
*  4.6.2 getting torsion parameters from FF file
	ITT=0
	do I=1,NT
	  IT=IT1(I)
	  JT=IT2(I)
	  KT=IT3(I)
	  LT=IT4(I)
	  J0=1
	  IFIND=0
* 4.6.2.1 Looking for exact dihedrals
	  ISUC=0
  	  do J=J0,NTOR
	    if((NFF(IT).eq.NT1(J).and.
     +NFF(JT).eq.NT2(J).and.NFF(KT).eq.NT3(J).and.NFF(LT).eq.NT4(J))
     +.or.(NFF(LT).eq.NT1(J).and.NFF(KT).eq.NT2(J).and.
     +NFF(JT).eq.NT3(J).and.NFF(IT).eq.NT4(J)))then
	       IFIND=IFIND+1
	       ITT=ITT+1
	       ITUK(ITT)=J
	       ITUKA(ITT)=I
	       ISUC=1
            end if
	  end do
*  4.6.2.2 Now looking for wildcards
	  if(ISUC.eq.0)then
  	    do J=J0,NTOR
	      if(((NFF(IT).eq.NT1(J).or.NT1(J)(1:1).eq.'X').and.
     +NFF(JT).eq.NT2(J).and.NFF(KT).eq.NT3(J).and.(NFF(LT).eq.NT4(J).
     +or.NT4(J)(1:1).eq.'X')).or.((NFF(LT).eq.NT1(J).or.NT1(J)(1:1)
     +.eq.'X').and.NFF(KT).eq.NT2(J).and.NFF(JT).eq.NT3(J).and.
     +(NFF(IT).eq.NT4(J).or.NT1(J)(1:1).eq.'X')))then
	        IFIND=IFIND+1
	        ITT=ITT+1
	        ITUK(ITT)=J
	        ITUKA(ITT)=I
              end if
	    end do
	  end if
	  if(IFIND.eq.0)then
	    ITT=ITT+1
	    ITUK(ITT)=0
	    ITUKA(ITT)=I
  	    write(*,*)' parameters for torsion ',NM(IT),'-',NM(JT),'-',
     +NM(KT),'-',NM(LT),' or ',NFF(IT),'-',NFF(JT),'-',NFF(KT),
     +'-',NFF(LT),' not found'
  	    write(20,*)' parameters for torsion ',NM(IT),'-',NM(JT),'-',
     +NM(KT),'-',NM(LT),' or ',NFF(IT),'-',NFF(JT),'-',NFF(KT),
     +'-',NFF(LT),' not found'
	  end if
	  if(IFIND.gt.1)then
  	    write(*,*)IFIND,' torsions ',NM(IT),'-',NM(JT),'-',
     +NM(KT),'-',NM(LT),' found'
  	    write(20,*)IFIND,' torsions ',NM(IT),'-',NM(JT),'-',
     +NM(KT),'-',NM(LT),' found'
	  end if
	  NTP=NTP+IFIND
	end do
*  4.6.3 Writing torsions to .mmol file
	write(13,'(a)')'#  Number of torsions'
	if(ISKIP.eq.1)then
	   write(13,*)NTP
	else
	   write(13,*)ITT
	end if
	do J=1,ITT
	   I=ITUKA(J)
	   IP=ITUK(J)
	   if(IP.gt.0)then
 	     write(13,'(4I6,2f9.3,I4,2x,4(a1,a4))')
     +   IT1(I),IT2(I),IT3(I),IT4(I),RT(IP),FT(IP),iabs(NMUL(IP)),
     +'#',NFF(IT1(I)),'-',NFF(IT2(I)),'-',NFF(IT3(I)),'-',NFF(IT4(I))
	   else
 	     if(ISKIP.ne.1)write(13,'(a1,4I6,a,2x,4(a1,a4))')
     +'!',IT1(I),IT2(I),IT3(I),IT4(I),
     +'   0.0     0.0    1 ',
     +'#',NFF(IT1(I)),'-',NFF(IT2(I)),'-',NFF(IT3(I)),'-',NFF(IT4(I))
	   end if
	end do
	write(*,*)' torsions list completed'
	write(20,*)' torsions list completed'
*  4.7 Impropers
C  TODO
*  4.8  Special 1-4 interactions
*  4.8.1  Finding pairs
	IP14=0
*  1-4 pairs from the list of torsions
	do I=1,NT
	   IT=IT1(I)
	   LT=IT4(I)
	   do M=1,NP14
	     if((NPI(M).eq.NFF(IT).and.NPJ(M).eq.NFF(LT)).or.
     &             (NPJ(M).eq.NFF(IT).and.NPI(M).eq.NFF(LT)))go to 191
	   end do
	   go to 195
*      pair found
 191	   IP14=IP14+1
	   IU14(IP14)=IT
	   JU14(IP14)=LT
	   MU14(IP14)=M
	   write(*,'(a,4i6,2a4)')
     + ' sp14 pair:',IP14,IT,LT,M,NM(I),NM(J)
 195	   continue
	end do 
	write(*,*)IP14,' special 1-4 pairs found'
	write(20,*)IP14,' special 1-4 pairs found'
	if(NSP14.gt.0)then
	  write(13,'(a)')'#   Special 1-4 interactions '
	  write(13,'(a)')'special'
	  write(13,*)NSP14
	  write(13,'(a)')'#  I   Sig14   Eps14'
	  do I=1,NSTS
	     J=IUK14(I)
	    if(J.gt.0)write(13,'(i6,2f11.4,a,a4)')
     + I,SGM14(J),EPS14(J),'   #  ',NFF(I)
	  end do
	  write(*,*)' list of ',NSP14,' special 1-4 parameters added'
	  write(20,*)' list of ',NSP14,' special 1-4 parameters added'
	end if
*  4.9  Special 1-4 pairs
	if(IP14.gt.0)then
	  write(13,'(a)')'#  Special 1-4 pairs'
	  write(13,'(a)')'par14'
	  write(13,*)IP14
	  write(13,'(a)')'#   I    J     Sig14    Eps14'
	  do K=1,IP14
	     I=IU14(K)
	     J=JU14(K)
	     M=MU14(K)
	     write(13,'(2i7,2f10.4)')I,J,SGP14(M),EPP14(M)
	  end do
	  write(*,*)' list of ',IP14,' special 1-4 pairs added'
	  write(20,*)' list of ',IP14,' special 1-4 pairs added'
	end if
*  4.10  1-4 scaling factors
	if(fudgeQQ.ne.1.)write(13,'(a,f8.4)')'sel14 ',fudgeQQ
*  5. Finish
	write(*,*)' Impropers, if needed,  should be added manually'
	close(13)
 	stop
*  ...
  99	IE=99
	STR=TAKESTR(12,IE)
  199	IE=99
	STR=TAKESTR(11,IE)
        write(*,*)' (..|..) '
 299	write(*,*)' Error in the list of bonds. Bond ',I
	end
*
*================ LENG ===============================================
*
	function LENG(STR,LS)
	character*1 STR(LS)
	IL=0 
	do I=1,LS
	  ISTR=ichar(STR(I))
	  if(ISTR.le.15)then
	    LENG=I-1
	    return
	  end if
        if(STR(I).ne.' ')IL=I
	end do 
	LENG=IL
	return
      end    
C
C============== TAKESTR =========================================
C
*   3.1 Definitions
      character*128 function TAKESTR(KAN,IE)
      character*128 AUX, fn 
      integer LCOUNT(256)   
      data LCOUNT /256*0/
      save AUX,LCOUNT
*   3.2 Normal action           
      if(IE.eq.99)go to 10
C  LCOUNT variable count line number (including commentaries) 
 1	LCOUNT(KAN)=LCOUNT(KAN)+1
 	read(KAN,'(a128)',err=10,end=20)AUX
	if(AUX(1:1).eq.'#'.or.AUX(1:1).eq.'*')go to 1
	TAKESTR=AUX 
	IE=0
	return
*   3.3 Error diagnostic
 10	write(*,*)'!!! error in input file '   
	if(KAN.eq.5)then
	  write(*,*)' Standard input ' 
	else
	  inquire(unit=KAN,name=fn) 
	  write(*,*)' File ',fn
	end if
	write(*,*)' in line ',LCOUNT(KAN)
	write(*,*)AUX
	stop
*    3.4 End of file case (is not signaled if input IE < 0 )
 20	if(IE.ge.0)then
 	  write(*,*)'!!! end of file reached' 
	  if(KAN.eq.5)then
	    write(*,*)' Standard input ' 
	  else
	    inquire(unit=KAN,name=fn) 
	    write(*,*)' File ',fn
	  end if
	  stop
	end if  
	TAKESTR='  ' 
	IE=-1
	return
	end 
