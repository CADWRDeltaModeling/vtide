

      module constituent

      integer :: ntotal = 0  ! total # constituents considered, including shallow water
      integer :: ntidal
      parameter(MC=70,mc2=mc*2,nmaxall = 170,nsatellite = 180,nshallowparent=320)
      integer :: II(MC2),JJ(MC2),KK(MC2),LL(MC2),MM(MC2),NN(MC2)
      real*8  :: SEMI(MC2)  ! todo: double check should be real*8
      character*5 :: kontab(MC)
      real*8      :: freq(MC)
 
      character*5 :: konco(nshallowparent) !todo: correct data type?
      real*8      :: coef(nshallowparent)

      character*5, :: kon(nmaxall)
      real*8 V(nmaxall),U(nmaxall),F(nmaxall)
      integer    :: NJ(nmaxall)
      real*8  :: pi=4.d0*atan(1.d0)
      real*8 :: twopi=8.d0*atan(1.d0)
      
      
      real*8 EE(nsatellite),PH(nsatellite)
      integer :: LDEL(nsatellite),MDEL(nsatellite)
      integer :: NDEL(nsatellite),IR(nsatellite)
      
      contains
 


      SUBROUTINE VUF(konx,vx,ux,fx)

      implicit none
C***********************************************************************
C*  THIS SUBROUTINE CALCULATES V (ASTRONOMICAL PHASE ARGUMENT), U AND F
C*  (NODAL MODULATION PHASE AND AMPLITUDE CORRECTIONS) FOR ALL CONSTITU-
C*  ENTS.
c*
C*	This October 1992 version also recalculates the constituent
c	frequencies for the middle of the analysis period
C

      real*8 ux,vx,fx

      integer k
      character*5 KONX
      integer :: lp= 6 ! todo: many redundant entries
      
C
C***********************************************************************
C*  THE DIMENSION OF KON, VU, F, AND NJ SHOULD BE AT LEAST EQUAL TO THE
C*  TOTAL NUMBER OF POSSIBLE CONSTITUENTS (PRESENTLY 146), THE DIMENSION
C*  OF II, JJ, KK, LL, MM, NN AND SEMI SHOULD BE AT LEAST EQUAL TO THE
C*  NUMBER OF MAIN CONSTITUENTS (PRESENTLY 45), THE DIMENSION OF EE,
C*  LDEL, MDEL, NDEL, IR, AND PH SHOULD BE AT LEAST EQUAL TO THE TOTAL
C*  NUMBER OF SATELLITES TO ALL THE MAIN CONSTITUENTS PLUS THE NUMBER
C*  OF CONSTITUENTS WITH NO SATELLITES (PRESENTLY 162+8),
C*  AND THE DIMENSION OF KONCO, AND COEFF SHOULD BE AT LEAST EQUAL TO
C*  THE SUM FOR ALL SHALLOW WATER CONSTITUENTS OF THE NUMBER OF MAIN
C*  CONSTITUENTS FROM WHICH EACH IS DERIVED (PRESENTLY 251).
C***********************************************************************
C* GIVEN CONSTITUENT KONX , THE NODAL CORRECTIONS V+U AND F ARE RETURNED
C
      DO 20 K=1,NTOTAL
      IF(KON(K).eq.KONX) go to 40
20    CONTINUE
      WRITE(LP,30)KONX
30    FORMAT('VUF STOP.',A5)
      STOP
40    VX=V(K)
      ux=u(k)
      FX=F(K)
      RETURN
      end subroutine
C
C***********************************************************************
C*  THE ASTRONOMICAL ARGUMENTS AND THEIR RATES OF CHANGE,
C*  S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP,  ARE READ FROM TWO RECORDS IN
C*  THE FORMAT(5F13.10):
C*     S0  = MEAN LONGITUDE OF THE MOON (CYCLES) AT 000 ET 1/1/1976.
C*     H0  = MEAN LONGITUDE OF THE SUN.
C*     P0  = MEAN LONGITUDE OF THE LUNAR PERIGEE.
C*     ENP0= NEGATIVE OF THE MEAN LONGITUDE OF THE ASCENDING NODE.
C*     PP0 = MEAN LONGITUDE OF THE SOLAR PERIGEE (PERIHELION).
C*     DS,DH,DP,DNP,DPP ARE THEIR RESPECTIVE RATES OF CHANGE OVER A 365
C*     DAY PERIOD AS OF 000 ET 1/1/1976.
C
      subroutine vt_read_constituent_db(name_constituent_db)

      implicit none
      character*80 :: name_constituent_db
      integer, parameter :: MC=70     !todo: many redundant entries
      integer, parameter :: MC2=MC*2
      integer :: kr=8
      
      real*8 dpp,ds,dp,dnp,dh
      !real*8 p, tau, dtau
      real*8 s0,h0,p0,enp0,pp0
      
      integer :: jbase,j1,jl ! todo: more global than this?
      integer :: j,k,k1,j4
      !integer index(nmaxall)
      open(unit=KR,file=name_constituent_db)

c	These values are no longer used though they are still
c	read in. More accurate polynomial approximations are 
c	now employed.
      READ(KR,50)S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP
 50   FORMAT(5F13.10)
C
C***********************************************************************
C*  HERE THE MAIN CONSTITUENTS AND THEIR DOODSON NUMBERS ARE READ IN
C*  FORMAT (6X,A5,1X,6I3,F5.2,I4). THE VALUES ARE RESPECTIVELY
C*     KON    = CONSTITUENT NAME
C*  II,JJ,KK,LL,MM,NN = THE SIX DOODSON NUMBERS
C*     SEMI   = PHASE CORRECTION
C*     NJ     = THE NUMBER OF SATELLITES FOR THIS CONSTITUENT.
C*  THE END OF ALL MAIN CONSTITUENTS IS DENOTED BY A BLANK CARD.
C
      JBASE=0
      DO 90 K=1,1000
      READ(KR,60)KON(K),II(K),JJ(K),KK(K),LL(K),MM(K),NN(K),SEMI(K),NJ(K)
      IF(len_trim(KON(K)) == 0) go to 100
70    J1=JBASE+1
      IF(NJ(K).GE.1) GO TO 75
      NJ(K)=1
      JL=J1
      PH(J1)=0.
      EE(J1)=0.
      LDEL(J1)=0
      MDEL(J1)=0
      NDEL(J1)=0
      IR(J1)=0
      GO TO 90
75    JL=JBASE+NJ(K)
C
C***********************************************************************
C*  IF NJ>0, INFORMATION ON THE SATELLITE CONSTITUENTS IS READ , THREE
C*  SATELLITES PER CARD, IN THE FORMAT (11X,3(3I3,F4.2,F7.4,1X,I1,1X)).
C*  FOR EACH SATELLITE THE VALUES READ ARE
C*     LDEL,MDEL,NDEL = THE CHANGES IN THE LAST THREE DOODSON NUMBERS
C*                  FROM THOSE OF THE MAIN CONSTITUENT.
C*     PH     = THE PHASE CORRECTION
C*     EE     = THE AMPLITUDE RATIO OF THE SATELLITE TIDAL POTENTIAL TO
C*            THAT OF THE MAIN CONSTITUENT.
C*     IR     = 1 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
C*            LATITUDE CORRECTION FACTOR FOR DIURNAL CONSTITUENTS
C*            2 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
C*            LATITUDE CORRECTION FACTOR FOR SEMI-DIURNAL CONSTI-
C*            TUENTS.
C*            OTHERWISE IF NO CORRECTION IS REQUIRED TO THE AMPLITUDE
C*            RATIO.
C
      READ(KR,80)(LDEL(J),MDEL(J),NDEL(J),PH(J),EE(J),IR(J),J=J1,JL)
90    JBASE=JL
60    FORMAT(6X,A5,1X,6I3,F5.2,I4)
80    FORMAT((11X,3(3I3,F4.2,F7.4,1X,I1,1X)))
100   NTIDAL=K-1
C
C***********************************************************************
C*  THE SHALLOW WATER CONSTITUENTS AND THE MAIN CONSTITUENTS FROM WHICH
C*  THEY ARE DERIVED ARE READ IN HERE WITH THE FORMAT
C*  (6X,A5,I1,2X,4(F5.2,A5,5X)). THE VALUES ARE RESPECTIVELY
C*     KON    = NAME  OF THE SHALLOW WATER CONSTITUENT
C*     NJ     = NUMBER OF MAIN CONSTITUENTS FROM WHICH IT IS DERIVED.
C*     COEF,KONCO = COMBINATION NUMBER AND NAME OF THESE MAIN
C*              CONSTITUENTS.
C*  THE END OF THESE CONSTITUENTS IS DENOTED BY A BLANK CARD.
C
      JBASE=0
      K1=NTIDAL+1
      DO 160 K=K1,1000
      J1=JBASE+1
      J4=J1+3
      READ(KR,130)KON(K),NJ(K),(COEF(J),KONCO(J),J=J1,J4)
      IF(len_trim(KON(K)) == 0) go to 170
160   JBASE=JBASE+NJ(K)
130   FORMAT(6X,A5,I1,2X,4(F5.2,A5,5X))
170   NTOTAL=K-1
      RETURN
      end subroutine
C
C***********************************************************************
C*  NTIDAL IS THE NUMBER OF MAIN CONSTITUENTS
C*  NTOTAL IS THE NUMBER OF CONSTITUENTS (MAIN + SHALLOW WATER)
C*  FOR  THE GIVEN TIME hr, THE TABLE OF F AND V+U VALUES IS
C*  CALCULATED FOR ALL THE CONSTITUENTS.
C*     F IS THE NODAL MODULATION ADJUSTMENT FACTOR FOR AMPLITUDE
C*     U IS THE NODAL MODULATION ADJUSTMENT FACTOR FOR PHASE
C*     V IS THE ASTRONOMICAL ARGUMENT ADJUSTMENT FOR PHASE.
c
c	setvuf calculates the V,u,f values at time hr for all constituents
C
      subroutine SETVUF(hr,xlat,mf,name)
      !use analysis_request ! todo this dependency should be changed to an argument
      implicit none
      integer :: mf
      character*5, dimension(mf) :: name
      !integer mc,mc2
      !parameter (MC=70,MC2=MC*2)
      real*8 :: xlat
      real*8 d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp,hh,tau,dtau,hr

      integer indx(nmaxall)
      integer kd0,int24,intdys,jbase,k,L
      integer iv,j1,jL,iuu
      integer j
      real*8 vdbl ! todo: Correct type?
      real*8 slat,sumc,sums,rr,uudbl,uu,vv
      integer k1,iflag,lk,km1,lp
     
 
      SLAT=SIN(PI*XLAT/180.)
c      CALL GDAY(1,1,76,19,KD)
c      YEARS=(hr/24.D0-KD)/365.00D0
C
C***********************************************************************
C*  THE ASTRONOMICAL ARGUMENTS ARE CALCULATED BY LINEAR APPROXIMATION
C*  AT THE MID POINT OF THE ANALYSIS PERIOD.
C
c      S=S0+YEARS*DS
c      H=H0+YEARS*DH
c      P=P0+YEARS*DP
c      ENP=ENP0+YEARS*DNP
c      PP=PP0+YEARS*DPP
c	day number measured from January 0.5 1900 (i.e.,
c	1200 UT December 31, 1899
      d1=hr/24.d0
      call gday(31,12,99,18,kd0)
      d1=d1-dfloat(kd0)-0.5d0
      call astr(d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp)
      INT24=24
      INTDYS=int((hr+0.00001)/INT24)
      HH=hr-dfloat(INTDYS*INT24)
      TAU=HH/24.D0+H-S
      dtau=365.d0+dh-ds
C
C***********************************************************************
C*  ONLY THE FRACTIONAL PART OF A SOLAR DAY NEED BE RETAINED FOR COMPU-
C*  TING THE LUNAR TIME TAU.
C
      JBASE=0
      DO 210 K=1,NTIDAL
      do 209 l=1,mf    !todo make sure we have access
      if(kon(k).eq.name(l)) then
c      FREQ(l)=(II(K)*DTAU+JJ(K)*DS+KK(K)*DH+LL(K)*DP+MM(K)*DNP+
c     1NN(K)*DPP)/(365.*24.)
      indx(k)=l
      end if
209   continue
      VDBL=II(K)*TAU+JJ(K)*S+KK(K)*H+LL(K)*P+MM(K)*ENP+NN(K)*PP+SEMI(K)
      IV=VDBL
      IV=(IV/2)*2
      Vv=VDBL-IV            ! rounded to 2?
      J1=JBASE+1
      JL=JBASE+NJ(K)
      SUMC=1.
      SUMS=0.
      DO 200 J=J1,JL
C
C***********************************************************************
C*  HERE THE SATELLITE AMPLITUDE RATIO ADJUSTMENT FOR LATITUDE IS MADE
C
      RR=EE(J)
      L=IR(J)+1
      GO TO (901,902,903),L
  902 RR=EE(J)*0.36309*(1.-5.*SLAT*SLAT)/SLAT
      GO TO 901
  903 RR=EE(J)*2.59808*SLAT
  901 CONTINUE
      UUDBL=LDEL(J)*P+MDEL(J)*ENP+NDEL(J)*PP+PH(J)
      IUU=UUDBL
      UU=UUDBL-IUU
      SUMC=SUMC+RR*COS(UU*TWOPI)
      SUMS=SUMS+RR*SIN(UU*TWOPI)
  200 CONTINUE
      F(K)=SQRT(SUMC*SUMC+SUMS*SUMS)
      v(k)=vv
      U(K)=ATAN2(SUMS,SUMC)/TWOPI
210   JBASE=JL
C
C***********************************************************************
C*  HERE F AND V+U OF THE SHALLOW WATER CONSTITUENTS ARE COMPUTED FROM
C*  THE VALUES OF THE MAIN CONSTITUENT FROM WHICH THEY ARE DERIVED.
C
      JBASE=0
      K1=NTIDAL+1
      IF(K1.GT.NTOTAL) RETURN
      DO 270 K=K1,NTOTAL
      F(K)=1.0
      V(K)=0.0
      u(k)=0.
      iflag=0
      do 269 lk=1,mf
      if(kon(k).eq.name(lk)) then
c      FREQ(lk)=0.
      iflag=1
      go to 268
      end if
269   continue
268   J1=JBASE+1
      JL=JBASE+NJ(K)
      DO 260 J=J1,JL
      KM1=K-1
      DO 240 L=1,KM1
      IF(KON(L).eq.KONCO(J)) go to 250
240   CONTINUE
      WRITE(LP,241)KONCO(J)
241   FORMAT('  SETVUF STOP.',A5)
      STOP
250   F(K)=F(K)*F(L)**ABS(COEF(J))
      V(K)=V(K)+COEF(J)*V(L)
      U(K)=U(K)+COEF(J)*U(L)
c      if(iflag.eq.1) FREQ(lk)=FREQ(lk)+COEF(J)*FREQ(indx(L))
260   continue
270   JBASE=JL
      RETURN
      END subroutine




      end module
