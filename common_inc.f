C  MDynaMix v.5.3
C
C  This is the include file containing static dimensioning
C  of the working fields. Practically all modules depend somehow on it, 
C  so you must recompile the sources each time you change this file.
C
      IMPLICIT REAL*8 (A-H,O-Z) 
*
*   1. Define maximum sizes of working arrays
*   -----------------------------------------
C   This part is now in a separate file dimpar.h
      include "dimpar.h"
*
*
*   2. Constants
*   ------------
C   Specfied here are constants available in all program modules which
C   have this file included
C
C   Physical constants:
C   Avagadro number, Boltzmann constant, J/cal conversion factor:
      PARAMETER (AVSNO=6.02252D23,BOLTZ=1.38054D-23,FKCALJ=4.1868)
C   Electron charge, Vaccuum "dielectric" constant:  
      PARAMETER (ELECHG=1.602D-19,EPS0=8.854D-12)
C   Planck constant
      PARAMETER (HPLANK=1.0545d-34)
C   Parameters to calculate erfc function      
      PARAMETER (A1=0.254829592D0,A2=-0.284496736D0,A3=1.421413741D0)
      PARAMETER (A4=-1.453152027D0,A5=1.061405429D0,B1=0.3275911D0)
C   Pi and angles      
      PARAMETER (PI=3.1415926536D0,TOTPI=1.128379169D0,PI2=2.D0*PI)
      PARAMETER (TODGR=180.D0/PI,ODIN=1.d0)
      parameter (cos30=0.8660254037844386,BOXYC=0.5/cos30,SQ3=2.*cos30)
*      
*   3. Variables
*   ------------
C   Specfied here are common variables available in all program modules which
C   have this file included
*
*   3.1. Names (character)
C
C   NM     - name of sites
C   NAMOL  - name of molecules;  NAME - short molecule name
C   PATHDB - path to molecular database
C   Lenght of NAMOL should be sinchronized with p.5.2, input.f
      CHARACTER NM*4,NAMOL*64,PATHDB*112,FFT*8
      CHARACTER*8 NAME
C   label -  program version
C   File names:
C        fdump  - restart file
C        ftrk   - Trajectory file
C        fmol   - Coordinate file for molecular viewers
C        filref - Some atoms may be bind to some locations, specified in 
C                 this file 
      character*32 label,fdump,frdf,fpairs,ftrk,fmol,filref,fname
      common /NAMES/label,fdump,ftrk,fname,fpairs,fmol,
     + PATHDB,filref
C       masses: MASS - at. units; MASSD - dimensional; MASSDI - inverse
*
*    3.2 Variables
*
C    More information on used variables are given in appropriate parts of program
C
      REAL*8 MASS,MASSD,MASSDI
C   Integers  
C   Some important numbers:
C         NOP    - total molecules
C         NSITES - total sites
C         NTYPES - total molecule types
C         other - see description in appropriate part of the source
      COMMON /INTEG/ NOP,NSITES,NTYPES,NSTOT,NFREQ,MOLINT,IPRINN
     X               ,NSTEP,NSTEPS,MSTEP,NSTTOT,NNAV,NAVT,NHIST,IHIST
     X               ,NAVER,IPRINT,NTRAJ,NSTEPA,NDUMP,IAVER,ICHECK
     +               ,ITREK,NTREK,NTRKF,LTREK,NTRKZ,ICFN,ICHNB,KMAX  
     +	     ,NRBT,NRAT,NRRT,NRIT,NRAA,NRBB,NRTT,NRII,NHDM,IO,NKV,
     +        NNAD,NRBND,NTYP1,ICELL,NPULL,IEXT,ISTAV,ICMOM,
     +        NNBT,NBEADS,JBEAD,NEKS,IPIMA,JNH,ICR,IVSLOW,ITSLOW,
     +        IPASS,NPASS,MPASS,IWLI,NWLI,NNADP,IOUTPUT,IDEE
C   energy/temperature related variables     
      COMMON /ENERGY/ PE,TEMPF,ENERF,COULF,UNITE,UNITM,UNITT
     X               ,UNITP,TOTMAS,DT,TSTEP,TSTEB,HTSTEP,HTSTEB,RHO
     X               ,FNST,TKE,TRTEMP,TDELT,TRPRES,FEXP,TOLER,BDTOL
     +               ,TRYCK,TRID,TRYCKM,PINT,CONVET,RTP,TKEQ,BETA
     X               ,TEMP,FTIMES,FTSCL,FNOP,RCUTT,EXAMPL,EXFREQ	
     +               ,EABSS,EABS,EFIELD,POTEN,DQQT,BTAU,QEK,QEKS,FNH
     +               ,TTR,TROT,TINT,FNSTI,FNSTR,RCUT,TEMPV,FACVT,RLR
     +               ,PE1,PELS1,PE2,PELS2,SPE,TTREK,EFACT,PERMOL,EPD	
     +		     ,TRYCKX,TRYCKY,TRYCKZ,RFF,RFF2,FBMASS,CONVEQ,
     +          BEADFAC,TEMPNOD,SPEI,TRPRESX,TRPRESY,TRPRESZ,TOKJ
C   Virials
      common /VIRRR/ VIRB,VIRA(3),VIRAN(3),WIRS,WIRSS,VIR1,VIR2,VIRD,
     +               VIRX,VIRY,VIRZ,VIRFX,VIRFY,VIRFZ,VIRPI,
     +               VIRKX,VIRKY,VIRKZ   
C   Control NVT/NPT ensembles
      COMMON /NVTNPT/ QT,QPR,DQT,SC,SCV,DQP,SCL,SCLX,SCLY,SCLZ,
     +                SCM(NTPS)
C   Size related quantities
      COMMON /SIZES/  UNITL,BOXL,BOYL,BOZL,HBOXL,HBOYL,HBOZL,ALPHA
     X               ,ALPHAD,RSSQ,SHORT,SHORT2,SHORTA,BOXY3,RDNA,VOL
     X               ,BOXLT,BOYLT,BOZLT,TEMPT,CHKTF
C   These numbers keep CPU time for different part of the program     
      common /TIMESZ/ TIM,TIMA,time0,timel,timeg,timef,timee,timen,
     +timeb,timea,timet,timest,timev,times  
C   Logical constants      
      LOGICAL        LSCLT,LSCFT,LNVT,LRST,LDMP,LCHL,LFORMI,LFORMO    
     X               ,LNPT,LMOVE,LRR,LSEP,LOCT,LHEX,LBONDL,LCMOM
     X               ,LINPUT,L0AVS,L0VS,L0CPU,LSHEJK,LCHP,LCHT,LEQ
     X               ,LXMOL,LMNBL,LGR,LGDMP,LGRST,LCHECK_TEMP        
     +               ,LMOL1,LMOL2,LCOMD,LGATHER,LSCTD,LVISUAL
     +               ,LPIMD,LPITH,LP_XY,LMDEE,LECH,LVSLOW,LTSLOW,LSCFL
     +               ,LDRAG,LDRAGT,LPAIRS,LPARFIL,LL2MC
      COMMON /LOGLOG/ LSCLT,LSCFT,LNVT,LRST,LDMP,LCHL,LFORMI,LFORMO
     X               ,LNPT,LXMOL,LRR,LSEP,LOCT,LHEX,LBONDL,LCMOM
     X               ,LINPUT,L0AVS,L0VS,L0CPU,LSHEJK,LCHP,LCHT,LEQ
     X               ,LMOVE(NTPS),LMNBL,LCHECK_TEMP        
     +               ,LMOL1,LMOL2,LCOMD,LGATHER,LSCTD,LVISUAL
     +               ,LPIMD,LPITH,LP_XY,LMDEE,LECH,LVSLOW,LTSLOW
     +               ,LSCFL,LDRAG,LDRAGT(NTPS),LPAIRS,LPARFIL,LL2MC
*
*   5. Arrays
*   ---------
C   Specfied here are common arrays available in all program modules which
C   have this file included. Dimensions of these arrays are defined above
C   in Sec.1
C
*     Molecules:
C   X,Y,Z              - molecular center of mass coordinates
C   PX,PY,PZ           - momentum of the molecule
C   QX,QY,QZ           - angular molecular momemtum
C   DPX,...            - dipole momentum
C   UPX,...            - unit vector fixed by two selected atoms
C   XPL,...            - average coordinates of "linked" atoms
        COMMON /COMCOM/  X(NPART), Y(NPART), Z(NPART)
     X               ,PX(NPART),PY(NPART),PZ(NPART)
        common /PULLL/ XPL(NPLM),YPL(NPLM),ZPL(NPLM)   
        common /EWALD/ RKX(NKVM),RKY(NKVM),RKZ(NKVM)
C     Atoms:
C   SX,SY,SZ           - atom coordinates
C   VX,VY,VZ           - atom momenta 
C   FX,FY,FZ           - forces (sometimes used as temporary arrays)
C   GX,GY,GZ           - slow forces
C   HX,HY,HZ           - fast forces
C   WX,WY,WZ           - local atom coordinates (relatively molecular center-of-mass)
C   OX,OY,OZ           - temporary array for atom coordinates or velocities
C   MASSDI             - (mass)**(-1)  (internal units)
C   Q                  - charge
      COMMON /SITSITX/ SX(NTOT) ,SY(NTOT) ,SZ(NTOT)
      COMMON /SITSITV/ VX(NTOT) ,VY(NTOT) ,VZ(NTOT)
      COMMON /SITSITF/ FX(NTOT) ,FY(NTOT) ,FZ(NTOT)
      COMMON /SITSITG/ GX(NTOT) ,GY(NTOT) ,GZ(NTOT)
      COMMON /SITSITH/ HX(NTOT) ,HY(NTOT) ,HZ(NTOT)
      COMMON /SITSITW/ WX(NTOT) ,WY(NTOT) ,WZ(NTOT)
      COMMON /SITSITO/  OX(NTOT) ,OY(NTOT) ,OZ(NTOT)
      COMMON /SITSITM/ MASSDI(NTOT),Q(NTOT)
      COMMON /SITSITT/ THERMO(NTOT,INHM),TEMBS(NTOT)
     +               ,THERMS(NTOT,INHM)
C   This arrays handle lists of neighbours    
	COMMON /LISTNB/MBSH,MBLN,NBS1(NBSMAX),NBS2(NBSMAX),
     +NNBB(NTOT),INBB(NTOT,NBDMAX)
      common /LIST1/NBL1(NBLMAX)
      common /LIST2/NBL2(NBLMAX)
C   Addresses and references
C   See file units.f for details
	common /ADDREFS/ IADDR(NTPS+1),ISADR(NTPS+1),ISADDR(NTPS+1),
     +                 ITM(NPART),ITYPE(NTOT),NNUM(NTOT),NSITE(NTOT),
     +                 MDX(NTPS,NTPS),ITS(NS),INPL(NPLM)
C   Molecule type characteristics   
C         NSPEC  - total molecules of this type
C         NSITS  - total atoms in a molecule of this type
C         ISHEJK - which molecules (1) are constrained
      COMMON /MOLTYP/ NSPEC(NTPS),NSITS(NTPS),IINIT(NTPS)              
     X              ,LIST(NTPS),ISHEJK(NTPS),NRCON(NTPS)
     X              ,NRB(NTPS),NRA(NTPS),NRT(NTPS),NRI(NTPS)
     X              ,IPOT(NTPS),NRTSC(NTPS),NAME(NTPS),NAMOL(NTPS)
C    some arrays for intermediate data     
      common /SOMMOL/ POTES(NTPP),POTLJ(NTPP),POTE1(NTPP),POTL1(NTPP)
     X               ,PELRC(NTPP),VRLRC(NTPP),SPEE(NTPS),SPEEI(NTPS)
     X               ,EKIN (NTPS+4),TEMPR(NTPS+4),SELFPE(NTPS)
     X               ,SUMMAS(NTPS),C14LJ(NTPS),C14EL(NTPS)
     X               ,TFACT(NTPS),QSUM(NTPS),PIRR(NS),PIRS(NS)
     +           ,PES14(NTPS),PSR14(NTPS),PES141(NTPS),PSR141(NTPS)
     +           ,FDRAG(NTPS)
C   Site characteristics
C        R  - initial local site coordinates from .mol file
C        NM - names of atom sites
      COMMON /SITTES/ SIGMA(NS),EPSIL(NS),CHARGE(NS),R(NS,3)
     +        ,EPAD(NNADM),SIGAD(NNADM),SIGT(NNBTM),EPST(NNBTM),
     +        SIGAT(NNBTM),EPSAT(NNBTM),A6LJ(NNBTM,NNBTM),
     +        B12LJ(NNBTM,NNBTM),A6_14(NNBTM,NNBTM),B12_14(NNBTM,NNBTM),
     +        INBT(NS),ILJ(NNADM),MASS(NS),MASSD(NS),NM(NS),FFT(NS),
     +        A6LJ0(NNBTM,NNBTM),B12LJ0(NNBTM,NNBTM),SIGAP(NNADM),
     +        EPAP(NNADM),ILJP(NNADM),JLJP(NNADM)     
C   Everything on bonds      
      COMMON /BONBON/ BB(NB),EB(NB),RB(NB),FB(NB),EBOND(NTPS),BR(NB)
     +	       ,BE(NB),DB(NB),RO(NB),IADB(NTPS),IB(NB),JB(NB),ID(NB)
     +	       ,IBI(NBO),JBJ(NBO),IBUK(NBO)
C   Everything on angles
      COMMON /ANGANG/ AA(NA),EA(NA),RA(NA),FA(NA),EANGL(NTPS),AR(NA)
     X               ,AE(NA),IADA(NTPS),IA(NA),JA(NA),KA(NA)
     +	       ,IAI(NAO),JAJ(NAO),KAK(NAO),IAUK(NAO)
C   Everything on torsions and impropers     
      COMMON /TORTOR/ TT (NT),ET (NT),ETORS(NTPS),FT1(NT),FT2(NT)
     X               ,FT3(NT),FT4(NT),FT5 (NT),TR(NT),TE(NT)
     X               ,RT (NT),FT (NT),NMUL(NT)
     X               ,IADT(NTPS),IT(NT),JT(NT),KT(NT),LT(NT)
     X               ,ITORS(NT),ITI(NTO),JTJ(NTO),KTK(NTO),LTL(NTO),
     +                ITUK(NTO)
C   Rotation matrix for each molecule     
      COMMON /ROTROT/ RMX(NPART,3,3)
C   Accumulators for averages and history of averages     
      COMMON /AVRAWR/ AV(NRQS),AW(NRQS),HIST(NRQS,LHIST)
C   RDF-calculations related arrays/variables      
*X      COMMON /RDFRDF/RDFINC,DRDFI,RDFCUT,RDFCUT2,LGR,LGDMP,LGRST
*X     X               ,MAXRDF,MAX,IGR(MAXGR),JGR(MAXGR),IIN(MAXGR)
*X     X               ,NRDF,IRDF(MAXS,MAXGR),NMGR(0:MAXGR)
*X     +		         ,NGRI(NS),IGRI(NS,MXGRS),MGRI(NS,MXGRS)
*XC   TCF-related arrays/variables
*X      CHARACTER*6 NAMECF
*X      PARAMETER (MCF=MAXCF*NTPS) 
*X	real*4 CC1,CC2,CC3,CC4,CC5,CC6
*X      COMMON /TCFTCF/ CC1(MAXTCF),CFV(MCF)
*X     X               ,CC2(MAXTCF),CFA(MCF)
*X     X               ,CC3(MAXTCF),CP1(MCF),CP2(MCF)
*X     X               ,CC4(MAXTCF),CU1(MCF),CU2(MCF)
*X     +               ,CC5(MAXTCF),CC6(MAXTCF)
*X     X               ,CVX(MCF),CVY(MCF),CVZ(MCF)
*X     X               ,CAX(MCF),CAY(MCF),CAZ(MCF)
*X     X               ,CCF(MCF)
*X     X               ,LCF,LCFDMP,LCFRST
*X     X               ,INX(MAXCF),NCALLS,NSTEG,NSTCF,NOM,JUMP
*X     X               ,N1(NTPS),N2(NTPS),ITCF(NTCF+6)
*X     X               ,NAMECF(NTCF)
*
*  6. Parameters for parallel execution (MPI)
*  ------------------------------------------
      PARAMETER (NBUFF  = 3*NTOT+2*NA+2*NB+2*NT)        
      integer NUMTASK,TASKID 
	integer NAB(0:NODES),NAE(0:NODES),NABS(0:NODES),NAP(0:NODES),
     +NAP3(0:NODES),NABS3(0:NODES)
      common /PARA/ NUMTASK,TASKID,IDG,LENMP,LENTT,LENIT,NBUF,NAP3,
     +NAP,NAB,NAE,NABS,NORD(NTPS),LENMB,LENMA,LENMT,LENMI,MAST,NAPC,
     +NABS3
*  Expanded ensemble
      logical LAUTOFACT
      common /MDEE/ EC(NENS),EE(NENS),EELRC(NENS),
     + AVE(NRQS,NENS),AWE(NRQS,NENS),PLOW,ENNO,WL_FAC,WL_INI,BEE,
     + WL_CORR,PW1,PW2,RAEX,LAUTOFACT
      common/MDEEI/ IWALK(NENS,NENS),IEE(NTPS),ITE(NENS,LHIST),
     +NE,ME,ME0,NDEL,NCEE,NNH,ICHE,ICHU
*    keep local variables in routines 
	save




