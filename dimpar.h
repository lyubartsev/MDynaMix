*   MDynaMix v 5.2
*
*   file dimpar.h
*
C     This file defines sizes of working arrays 
C     Edit parameters of this part to adjust to you SYSTEMS:
C     MOLECULAR and COMPUTER
C     Of course, the first one prefers as large as possible values and
C     the second one would like as low as possible values, so you should
C     compromise...       
C
C     Parameters which affect the memory most:
C     NPART, NTOT, MAXCF, NRQS, LHIST, NBLMX
C     (NBLMX scales as 1/NPROCS for parallel jobs)  
C
C     You must recompile the sources each time you change this file.
C
C
C  These parameters set maximum number of:              
C       NPART -  molecules
      PARAMETER (NPART  = 50000  )
C       NTPS - different types of molecules
      PARAMETER (NTPS   =   8  )
C       NS - sites (i.e., sum of atoms if one take one molecule of each type)
C       NTOT - total atoms
      PARAMETER (NS     = 6000  , NTOT   = 150000             ) 
C       NNBTM - maximum number of different non-bonded types
      PARAMETER(NNBTM=1000)
C       NBSMX - max neighbours per atom inside cutoff of fast forces (SHORT)
C       (not essential if SHAKE algorithm used)
C       NBLMX - max neighbours per atom inside cutoff Rcut
C       These parameters may be decreased for parallel jobs as 1/NPROCS
      PARAMETER (NBSMX=64, NBLMX=800)
      PARAMETER (NBSMAX = NBSMX*NTOT,NBLMAX=NBLMX*NTOT          )  
C       NB     different covalent bonds 
C             (goes over ALL bonds on ONE molecule of EACH type)  
C       NBO    total covalent bonds per processor
      PARAMETER (NB     = 20000  , NBO     = 200000             )
C       NKVM - maximum number of reciprocal vectors in the Ewald sum
      PARAMETER (NKVM = 100000)
C       NA,NAO  the same for covalent angles
      PARAMETER (NA    = 50000 , NAO     = 200000           )    
C       NT,NTO  the same for dihedral angles
      PARAMETER (NT     = 50000 , NTO   = 200000       )  
C       NNADM - number of specific Lennard-Jones parameters for 1-4
C               interactions 
C       NPLM  - max number of "linked" atoms 
C               (attached to some points by harmonic forces) 
      PARAMETER (NNADM=2000, NPLM=100 )                 
C              RDF resolution    num of RDFs     Max RDF for 1 site 
*X      PARAMETER (MAXS   = 1  , MAXGR  = 1 ,   MXGRS=1)        
C                    points for tcf calculations  ; total tcf
*XC   set MAXCF=1 if you do not want calculate tcf and want to save memory
*X      PARAMETER (MAXCF  = 1    , NTCF   = 12               )
*X      PARAMETER (MAXTCF = (MAXCF+2)*NPART*3                ) 
C                   different averages
      PARAMETER (NRQS   = 50000 , NTPP =(NTPS*(NTPS+1))/2+1 ) 
C   max number of "bound" neighbours 
      PARAMETER (NBDMAX    = 100 )           
C       intermediate averaging points
      PARAMETER (LHIST  = 500)                        
C       maximum number of nodes   
      parameter (nodes=512)
C       maximum number of ensembles in MDEE
      parameter (NENS=128)
C       maximum number of extra potentials
      parameter (MAXPOT=10,INHM=1)
*

