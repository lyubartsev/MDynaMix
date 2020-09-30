C    TRANAL v 5.0
C    Computation of average dielectric permittivity
C    from the fluctuation of the total dipole moment
C
      program DIEL
      include "tranal.h"
      call SETTRAJ
      read(*,*)TEMPR
      IEND=0
      DSX=0.
      DSY=0.
      DSZ=0.
      DS4=0.
      DS2=0.
      VOLS=0.
      do while(IEND.eq.0)
        call READCONF(IEND)
        call GETCOM
        DX=0.
        DY=0.
        DZ=0.
        ISP=0.
        IM=0
        do ITYP=1,NTYPES
          do J=1,NSPEC(ITYP)
            DMX = 0.
            DMY = 0.
            DMZ = 0.
            IM=IM+1
            do K=1,NSITS(ITYP)
               ISP=ISP+1
               IS=K+ISADR(ITYP)
               DMX=DMX+CHARGE(IS)*(SX(ISP)-X(IM))
               DMY=DMY+CHARGE(IS)*(SY(ISP)-Y(IM))
               DMZ=DMZ+CHARGE(IS)*(SZ(ISP)-Z(IM))
*               DMX=DMX+CHARGE(IS)*SX(ISP)
*               DMY=DMY+CHARGE(IS)*SY(ISP)
*               DMZ=DMZ+CHARGE(IS)*SZ(ISP)
            end do
            DX = DX + DMX
            DY = DY + DMY
            DZ = DZ + DMZ
          end do
        end do
        DSX=DSX+DX
        DSY=DSY+DY
        DSZ=DSZ+DZ
        D2=DX**2+DY**2+DZ**2
        DEB=sqrt(D2)
        DS2=DS2+D2
        DS4=DS4+D2**2
        VOLS=VOLS+VOL
        if(IPRINT.ge.6)write(*,'(a,f12.3,a,3f10.3)')
     +  ' time',FULTIM,' Dtot=',DEB,D2
      end do
      DSX=DSX/IAN
      DSY=DSY/IAN
      DSZ=DSZ/IAN
      D2S=DSX**2+DSY**2+DSZ**2
      DS2=DS2/IAN
      DS4=DS4/IAN
      DDD=3.*sqrt((DS4-DS2**2)/IAN)
      VOLS=VOLS/IAN
      write(*,*)IAN,' configurations analysed'
      write(*,*)' Average dipole moment:',DSX,DSY,DSZ
      write(*,*)'  <M>**2 = ',D2S
      write(*,*)'  <M**2> = ',DS2,' el**2*A**2   +/- ',DDD
      EPS = ELECHG**2*(DS2-D2S)/(EPS0*3.d-10*VOLS*BOLTZ*TEMPR)+1.
      DEPS = EPS*DDD/DS2
      write(*,*)' Volume ',VOLS,'   Cub ',VOLS**(1./3.)
      write(*,*)' Dielectric constant: ',EPS,'  +/- ',DEPS
      stop
      end
