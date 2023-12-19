C========================================
C
C     MDynaMix v.5.2
C
*     PART 11
*
*     File util.f
*     -------------
*
*============ GAUSS ===============================================
*
      real*8 FUNCTION GAUSS()
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (A1=3.949846138,A3=0.252408784)
      PARAMETER (A5=0.076542912,A7=0.008355968)
      PARAMETER (A9=0.029899776)
      SUM   = 0.0
      DO I  = 1,12
      SUM   = SUM + RANF()
      END DO                   
      R     = (SUM-6.0)/4.0
      R2    = R*R
      GAUSS = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R
      RETURN
      END
*
*============================== RANF =========================
*
      real*8 FUNCTION RANF()
      IMPLICIT real*8 (A-H,O-Z)
      INTEGER L,C,M,SEED
      PARAMETER (L=1029,C=221591,M=1048576)
      SAVE SEED
      DATA SEED /0/
      SEED    = MOD(SEED*L+C,M)
      RANF    = DFLOAT(SEED)/M
      RETURN	
      END               
*
*=============== RAND =======================================
*                             
*  Random number generator
*       subroutine RANDU(IX,IY,W)
        real*8 function RAND()
        implicit real*8 (A-H,O-Z)
        data ITWO/2/,M2/0/,IX/14371/
        save
        if(M2.ne.0)go to 20
          M=1
 10       M2=M
          M=ITWO*M2
          if(M.gt.M2)go to 10
          HALFM=M2
          IA0=HALFM*datan(1.d0)/8.d0
          IC0=HALFM*(0.5d0-dsqrt(3.d0)/6.d0)
          IA=8*IA0+5
          IC=2*IC0+1
          MIC=(M2-IC)+M2
          S=sngl(0.5d0/HALFM)
 20       IY=IX*IA
        if(IY.gt.MIC)IY=(IY-M2)-M2
        IY=IY+IC
        if(IY/2.gt.M2)IY=(IY-M2)-M2
        if(IY.lt.0)IY=(IY+M2)+M2
        IX=IY
        RAND=IY*S
        return
        end
*
*=============== RANQ =========================================
*
      SUBROUTINE RANQ(QP)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION QP(4)
*
      PI=3.14159D0
      A=PI*RANF()/2.D0
      B=PI*RANF()
      C=PI*RANF()
      QP(1)=DSIN(A)*DSIN(B-C)
      QP(2)=DSIN(A)*DCOS(B-C)
      QP(3)=DCOS(A)*DSIN(B+C)
      QP(4)=DCOS(A)*DCOS(B+C)
      RETURN
      END
*
*===================== INV3B3 ==========================================
*
      SUBROUTINE INV3B3(A,AI,DET,IER)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION A(9),AI(9)
*
      AI(1) = A(5)*A(9)-A(6)*A(8)
      AI(2) = A(3)*A(8)-A(2)*A(9)
      AI(3) = A(2)*A(6)-A(3)*A(5)
      AI(4) = A(6)*A(7)-A(4)*A(9)
      AI(5) = A(1)*A(9)-A(3)*A(7)
      AI(6) = A(3)*A(4)-A(1)*A(6)
      AI(7) = A(4)*A(8)-A(5)*A(7)
      AI(8) = A(2)*A(7)-A(1)*A(8)
      AI(9) = A(1)*A(5)-A(2)*A(4)
*
      DET   = A(1)*AI(1)+A(4)*AI(2)+A(7)*AI(3)
	SPUR=A(1)+A(5)+A(9)   
	IER=0
	if(dabs(DET).lt.1.d-6*dabs(SPUR**3))IER=1
      IF(DABS(DET).gt.0.D0) then
	  R=1.D0/DET
	else
	  R=0.d0
	  IER=2
	end if
*
      AI(1) = R*AI(1)
      AI(2) = R*AI(2)
      AI(3) = R*AI(3)
      AI(4) = R*AI(4)
      AI(5) = R*AI(5)
      AI(6) = R*AI(6)
      AI(7) = R*AI(7)
      AI(8) = R*AI(8)
      AI(9) = R*AI(9)
*
      RETURN
      END
*
*=============== ROTATE ==============================================
*
      SUBROUTINE ROTATE(Q,X,Y,Z)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION Q(4)
      SIGN = -1.D0
      AA   = Q(1)*Q(1)
      BB   = Q(2)*Q(2)
      CC   = Q(3)*Q(3)
      DD   = Q(4)*Q(4)
      AB   = Q(1)*Q(2)
      AC   = Q(1)*Q(3)
      AD   = Q(1)*Q(4)*SIGN
      BC   = Q(2)*Q(3)
      BD   = Q(2)*Q(4)*SIGN
      CD   = Q(3)*Q(4)*SIGN
      RX   = X*(-AA+BB-CC+DD)+Y*(+CD-AB)*2.D0 +Z*(+BC+AD)*2.D0
      RY   = X*(-AB-CD)*2.D0 +Y*(+AA-BB-CC+DD)+Z*(+BD-AC)*2.D0
      RZ   = X*(+BC-AD)*2.D0 +Y*(-AC-BD)*2.D0 +Z*(-AA-BB+CC+DD)
      X    = RX
      Y    = RY
      Z    = RZ
      RETURN
      END
*
*=============== GETMOI ========================================
*
      SUBROUTINE GETMOI(X,Y,Z,M,MOMENT,NSS)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 M,MOMENT
      DIMENSION X(NSS),Y(NSS),Z(NSS),M(NSS)
      DIMENSION MOMENT(3,3),R(3)
*
      DO INERT           = 1,3
      DO IA              = 1,3
      TT                 = 0.D0
      MOMENT(INERT,IA)   = 0.D0
      IF(INERT.EQ.IA) TT = 1.D0
      DO I               = 1,NSS
      R(1)               =  X(I)
      R(2)               =  Y(I)
      R(3)               =  Z(I)
      SQS                =  X(I)**2+Y(I)**2+Z(I)**2
      MOMENT(INERT,IA)   = MOMENT(INERT,IA)+M(I)*(TT*SQS-R(INERT)*R(IA))
      END DO! OF I
      END DO! OF IA
      END DO! OF INERT
*
      RETURN
      END
*
*========== HH3BY3 ============================================
*
      SUBROUTINE HH3BY3(A,D,E)
*
      IMPLICIT real*8 (A-H,O-Z)
*
      DIMENSION A(3,3),D(3),E(3)
      data ICOUNT/0/
*
      H       = 0.
      SCALE   = ABS(A(3,1))+ABS(A(3,2))
      IF(SCALE.EQ.0) THEN
      E(3)    = A(3,2)
      ELSE
      A(3,1)  = A(3,1)/SCALE
      A(3,2)  = A(3,2)/SCALE
      H       = A(3,1)**2+A(3,2)**2+H
      F       = A(3,2)
      G       =-SIGN(SQRT(H),F)
      E(3)    = SCALE*G
      H       = H-F*G
      A(3,2)  =   F-G
*
      A(1,3)  = A(3,1)/H
      G       = A(1,1)*A(3,1)+A(2,1)*A(3,2)
      E(1)    = G/H
      F       = E(1)*A(3,1)
*
      A(2,3)  = A(3,2)/H
      G       = A(2,1)*A(3,1)+A(2,2)*A(3,2)
      E(2)    = G/H
      F       = F+E(2)*A(3,2)
*
      HH      = F/(H+H)
      F       = A(3,1)
      G       = E(1)-HH*F
      E(1)    = G
      A(1,1)  = A(1,1)-F*E(1)-G*A(3,1)
      F       = A(3,2)
      G       = E(2)-HH*F
      E(2)    = G
      A(2,1)  = A(2,1)-F*E(1)-G*A(3,1)
      A(2,2)  = A(2,2)-F*E(2)-G*A(3,2)
      END IF
      D(3)    = H
*
      SCALE   = 0.
      E(2)    = A(2,1)
      D(2)    = 0.
*
      E(1)    = 0.
      D(1)    = A(1,1)
      A(1,1)  = 1.
      IF(D(2).NE.0) THEN
      G       = A(2,1)*A(1,1)
      A(1,1)  = A(1,1)-A(1,2)*G
      END IF
      D(2)    = A(2,2)
      A(2,2)  = 1.
      A(2,1)  = 0.
      A(1,2)  = 0.
      IF(D(3).NE.0) THEN
      G       = A(3,1)*A(1,1)+A(3,2)*A(2,1)
      A(1,1)  = A(1,1)-A(1,3)*G
      A(2,1)  = A(2,1)-A(2,3)*G
      G       = A(3,1)*A(1,2)+A(3,2)*A(2,2)
      A(1,2)  = A(1,2)-A(1,3)*G
      A(2,2)  = A(2,2)-A(2,3)*G
      END IF
      D(3)    = A(3,3)
      A(3,3)  = 1.
      A(3,1)  = 0.
      A(1,3)  = 0.
      A(3,2)  = 0.
      A(2,3)  = 0.
*
*Diagonalize it:
*
      E(1)    = E(2)
      E(2)    = E(3)
      E(3)    = 0.
      DO 15 L = 1,3
      ITER    = 0
    1 CONTINUE
      DO 12 M = L,2
      DD      = ABS(D(M))+ABS(D(M+1))
      IF(ABS(E(M))+DD.EQ.DD) GO TO 2
   12 CONTINUE
      M       = 3
    2 CONTINUE
      IF(M.NE.L) THEN
      IF(ITER.EQ.30) then
        write(*,*) '!!! TOO MANY ITERATIONS!    in HH3BY3'
        ICOUNT=ICOUNT+1
        if(ICOUNT.gt.1000)stop
        return
      end if
      ITER    = ITER+1
      G       = (D(L+1)-D(L))/(2.*E(L))
      R       = SQRT(G**2+1.)
      G       = D(M)-D(L)+E(L)/(G+SIGN(R,G))
      S       = 1.
      C       = 1.
      P       = 0.
      DO 14 I = M-1,L,-1
      F       = S*E(I)
      B       = C*E(I)
      IF(ABS(F).GE.ABS(G)) THEN
      C       = G/F
      R       = SQRT(C**2+1.)
      E(I+1)  = F*R
      S       = 1./R
      C       = C*S
      ELSE
      S       = F/G
      R       = SQRT(S**2+1.)
      E(I+1)  = G*R
      C       = 1./R
      S       = S*C
      END IF
      G       = D(I+1)-P
      R       =(D(I)-G)*S+2.*C*B
      P       = S*R
      D(I+1)  = G+P
      G       = C*R-B
      F       = A(1,I+1)
      A(1,I+1)= S*A(1,I)+C*F
      A(1,I)  = C*A(1,I)-S*F
      F       = A(2,I+1)
      A(2,I+1)= S*A(2,I)+C*F
      A(2,I)  = C*A(2,I)-S*F
      F       = A(3,I+1)
      A(3,I+1)= S*A(3,I)+C*F
      A(3,I)  = C*A(3,I)-S*F
   14 CONTINUE
      D(L)    = D(L)-P
      E(L)    = G
      E(M)    = 0.
      GO TO 1
      END IF
   15 CONTINUE
      RETURN
      END
*
*========== ORDER3x3=================================================
*
C   Set values of vector G in decreasing order and matrix A correspondingly
C
      subroutine ORDER3X3(A,G)
      real*8 A(3,3),G(3),B(3,3),H
      integer IO(3)
      IO(1)=1
      IO(2)=2
      IO(3)=3
      if(G(3).gt.G(2))then
         IO(2)=3                  !   G
         IO(3)=2
         H=G(3)
         G(3)=G(2)
         G(2)=H
         if(G(2).gt.G(1))then
           G(2)=G(1)
           G(1)=H
           IO(1)=3
           IO(2)=1
           if(G(3).gt.G(2))then
             H=G(3)
             G(3)=G(2)
             G(2)=H
             IO(2)=2
             IO(3)=1
           end if
         end if
       else
         if(G(2).gt.G(1))then
           H=G(1)
           G(1)=G(2)
           G(2)=H
           IO(1)=2
           IO(2)=1
           if(G(3).gt.G(2))then
             G(2)=G(3)
             G(3)=H
             IO(2)=3
             IO(3)=1
           end if
         end if
       end if
       if(G(1).lt.G(2).or.G(2).lt.G(3))write(*,*)' wrong G order:',G
       do I=1,3
         do J=1,3
           JJ=IO(J)
           B(I,J)=A(I,JJ)
         end do
       end do
       do I=1,3
          do J=1,3
             A(I,J)=B(I,J)
           end do
       end do
       return
       end
*....
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
        if(STR(I).ne.' '.and.ISTR.gt.16)IL=I
      end do 
      LENG=IL
      return
      end    
*                                                                              **================ DISP ============================================           *                                              
      SUBROUTINE DISP(FT,FM,FD,NI0,INN)
      implicit real*8 (A-H,O-Z)
      DIMENSION FT(INN)
*                                                                              *  Dispersion calculation                        
*                                                                      
      if(INN.le.NI0)then
        FM=FT(INN)
        FD=2.*FT(INN) 
        return  
      end if
      FM        = 0.0
      DO INU    = NI0,INN
        FM        = FM+FT(INU)
      END DO
      FM        = FM/(INN-NI0+1)
      FD        = 0.0
      DO INU    = NI0,INN
        FD        = FD+(FT(INU)-(FM))**2
      END DO
      FD        = DSQRT(FD)/(INN-NI0)
      RETURN
      END

