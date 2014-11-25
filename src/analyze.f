

!todo: no error trapping
!todo: make logging safe from other contexts instead of raw print

      subroutine vt_analyze(nobs,
     &                      itime,
     &                      obsdata,
     &                      latd,latm,  
     &                      mf,ndef,itrend, nloc,
     &                      name,freqc,
     &                      ninfer_cnstnts,ninfer,
     &                      infer_main_const,     
     &                      infer_names, 
     &                      infer_freq,     
     &                      infer_amp, 
     &                      infer_zeta,
     &                      ampc,phg,
     &                      ampci,phgi,
     &                      amaj,amin,ainc,g,
     &                      amaji,amini,ainci,gi,
     &                      eigval,
     &                      sig,
     &                      cov,
     &                      sdev,sdev0,
     &                      ssq,
     &                      rmsr,rmsr0,
     &                      resmax,resmax_ls,iresmax,
     &                      fitted,
     &                      istat)


      use constituent
      use file_io

      implicit none
      
      integer, parameter :: mcc = 70         ! max allowed consituent estimates
      integer, parameter :: nmaxp1 = mc*2
      integer, parameter :: max_infer = 80
      
      integer,intent(in) :: nobs             !< number of observations
      integer,intent(in) :: itime(7,nobs)    !< obs times ... year,mon,day,hour,min,sec,century  !todo: this is silly
      real*8,intent(in) :: obsdata(nobs,ndef,nloc) !< observations. Rows are times, columns are stations. For vectors station1_u, station1_v, station2_u, ...
      integer,intent(in) :: latd             !< representative latitude, degrees
      integer,intent(in) :: latm             !< representative longitude, minutes
      integer,intent(in) :: mf               !< number of main (ie, not inferred or satellite) frequencies in analysis
      integer,intent(in) :: ndef             !< number of componentes, 1 is scalar and 2 is for u,v pairs. other cases not handled 
      integer,intent(in) :: itrend           !< flag indicating a linear trend should be included: 0 = no, 1 = yes
      integer,intent(in) :: nloc
      character*5,intent(in),dimension(*) :: name    !< name of the main constituents in the analysis
      real*8,intent(in), dimension(*)     :: freqc   !< frequencies of the main constituents. Note these are not actually used because the astronomical argument is computed directly.
      integer,intent(in) :: ninfer_cnstnts    !< number of main constituents that have inferences attached to them
      integer,intent(in) :: ninfer(max_infer) !< number of inferences for each main constituent with inferences      
      character*5, dimension(max_infer),intent(in)   :: infer_main_const !< name of main constituent associated with each inference
      character*5,dimension(max_infer,16),intent(in) :: infer_names      !< name of inferred constituent associated with each inference
      real*8, dimension(max_infer,16) :: infer_freq           !< frequency of inferred constituent associated with each inference
      real*8, dimension(max_infer,16) :: infer_amp            !< amplitude ratio assumed for each inference
      real*8, dimension(max_infer,16) :: infer_zeta           !< zeta (phase relative to main constituent) for each inference

      
      real*8, intent(out) :: ampc(mcc,ndef,nloc)              !< amplitude estimate for each main contituent
      real*8, intent(out) :: phg(mcc,ndef,nloc)               !< phase estimate for each main constituent in degrees
      real*8, intent(out) :: ampci(10,10,ndef,nloc)           !< estimated amplitude of inferred constituents !todo: fix hardwired dimensions
      real*8, intent(out) :: phgi(10,10,ndef,nloc)            !< estimated phase of inferred constituents      
	real*8, intent(out) :: amaj(mcc)                        !< major ellipse parameter, main constituent            
	real*8, intent(out) :: amin(mcc)                        !< minor ellipse parameter, main constituent
	real*8, intent(out) :: ainc(mcc)                        !< incl angle, main constituent
	real*8, intent(out) :: g(mcc)                           !< phase, main constituent 
	real*8, intent(out) :: amaji(10,10)                     !< major ellipse param, inferred constituent
	real*8, intent(out) :: amini(10,10)                     !< minor ellipse param, inferred constituent
	real*8, intent(out) :: ainci(10,10)                     !< incl angle ellipse param, inferred constituent
	real*8, intent(out) :: gi(10,10)                        !< phase, inferred constituent
      real*8,intent (out) :: eigval(nmaxp1)                   !< eigenvalues     
      real*8,intent(out) ::  sig(4,nmaxp1,ndef,nloc)          !< significance statistics for each coefficient estimate, sig1= ... sig2 = ... sig3 = ... sig4=ttest
      real*8,intent(out) ::  cov(2*mf-1+itrend,2*mf-1+itrend)  !< covariance matrix of the variables in the least squares problem (z0 and linear trend if used plus 2 per constituent)

      real*8,intent(out)  ::  sdev(ndef,nloc)                   !< standard deviation !todo ... of what?
      real*8,intent(out)  ::  sdev0(ndef,nloc)                  !< standard deviation of the original right hand side
      real*8,intent(out) ::  ssq(ndef,nloc)                    !< sum squared error of each fit
      real*8,intent(out) ::  rmsr(ndef,nloc)                   !< root mean square error of each fit
      real*8,intent(out) ::  rmsr0(ndef,nloc)                  !< root mean square of ???
      real*8,intent(out) ::  resmax(ndef,nloc)                 !< maximume residual from each rhs
      real*8,intent(out) ::  resmax_ls(ndef,nloc)              !< maximume residual from each rhs as given by svd routine (should be same as resmax)
      integer,intent(out) :: iresmax(ndef,nloc)               !< index of the largest residual of each fit
      integer,intent(out) :: istat                            !< return code
      real*8,intent(out) :: fitted(nobs,ndef,nloc)
      
      real*8 :: amp(mcc)
      real*8 :: phi(mcc)

      real*8 :: toler
      real*8 :: wmax,wmin             ! used to calculate max and min eigenvalues, used to estimate conditioning of the problem
      real*8 av(ndef,nloc)            ! used for averages of rhs
      real*8 sum1,hr,hrm
      real*8 :: x(nobs)               ! time as a real*8 todo: used to be next to time
      real*8 :: time(nobs)            ! todo: wasteful variable ... couple scalars and x() are enough for time
      real*8 c,s
      real*8 sumtest
        ! whether to replace obs with fitted
c
C***********************************************************************
C*
C*  THIS PROGRAM DOES A TIDAL HEIGHTS 'HARMONIC' ANALYSIS OF IRREGULARLY
C*  SAMPLED OBSERVATIONS.  THE ANALYSIS METHOD IS A LEAST SQUARES FIT
C*  USING SVD COUPLED WITH NODAL MODULATION AND INFERENCE(IF SO REQUESTED).  
c*	
c*	The code is based on TOPEX analysis code originally developed by Josef 
c*	Cherniawsky (JAOT, 2001, 18(4): 649-664) and modified by Rob Bell and 
c*	Mike Foreman. Enhancements to that version include
c*
c*	1. Provision for multi-constituent inferences computed directly within
c*	  the least squares matrix rather than as post fit corrections. This
c*	  means that the inferred constituents will affect all constituents, not
c*	  just the reference constituent.
c*
c*	2. An extension to permit the analysis of current observations.
c*
c*	3. Removal of a central time as the basis for the calculation of thefi
c*	  astronomical arguments V. Now the V value for each observation, as well 
c*	  as those for the nodal corrections f and u (done for the JAOT analysis), 
c*	  are incorporated directly into the overdetermined matrix. These changes
c*	  mean that analyses no longer need be restricted to periods of a year or
c*	  less. (Though as the period approaches 18.6 years, using another "long 
c*	  period" analysis program that solves for the "nodal satellites" directly
c*	  is advisable.)
c* 
C***********************************************************************
C*
C*  FILE REFERENCE NUMBERS OF DEVICES REQUIRED BY THIS PROGRAM.
C*      KR  - INPUT FILE  - CONTAINS THE TIDAL CONSTITUENT INFORMATION.
C*      KR1 - INPUT FILE  - GIVES ANALYSIS TYPE AND TIDAL STATION
C*                          DETAILS.
C*      KR2 - INPUT FILE  - CONTAINS THE OBSERVED TIMES AND HEIGHTS.
C*      LP  - OUTPUT FILE - LINE PRINTER.
C*  PRESENTLY KR,KR1,KR2, AND LP ARE ASSIGNED THE RESPECTIVE VALUES
C*  8,5,9, AND 6.  SEE THE MANUAL OR COMMENT STATEMENTS WITHIN THIS
C*  PROGRAM FOR FURTHER DETAILS ON THEIR USE.
C*
C***********************************************************************
C*
C*  ARRAY DEFINITIONS AND DIMENSION GUIDELINES.
C*
C*  LET    MF      BE THE NUMBER OF MAIN CONSTITUENTS
C*         MCC      BE THE TOTAL NUMBER OF CONSTITUENTS, INCLUDING Z0
C*                 AND ANY INFERRED CONSTITUENTS, TO BE INCLUDED IN THE
C*                 ANALYSIS; (For T/P, MCC=30 > NUMBER OF CONSTITUENTS)
C*         NOBS    BE THE NUMBER OF OBSERVATIONS;
C*         NVAR    BE 2*MF-1;
C*         MEQ     BE NOBS*2 IF ALL THE OBSERVATIONS ARE EXTREMES AND
C*                 THE DERIVATIVE CONDITION IS TO BE INCLUDED FOR EACH,
C*                 AND NOBS OTHERWISE.
C*  THEN PARAMETERS NMAXP1, AND NMAXPM SHOULD BE AT LEAST NVAR+1, AND 
C*  NEQ+NVAR RESPECTIVELY.  THEY ARE CURRENTLY SET TO 40 AND 240.
C*
C*  NAME(I)        IS THE ARRAY CONTAINING ALL THE CONSTITUENT NAMES,
C*                 INCLUDING Z0 AND ANY INFERRED CONSTITUENTS, TO BE IN
C*                 THE ANALYSIS.  IT SHOULD BE DIMENSIONED AT LEAST MC.
C*  FREQC(I),      ARE THE ARRAYS OF FREQUENCIES IN CYCLES/HR AND
C*  FREQ(I)        RADIANS/HR RESPECTIVELY CORRESPONDING TO THE
C*                 CONSTITUENT NAME(I).  THEY SHOULD BE DIMENSIONED AT
C*                 LEAST MC.
C*  AMP(I),PHI(I)   ARE ARRAYS CONTAINING THE RAW AMPLITUDE AND PHASE FOR
C*                 CONSTITUENT NAME(I) AS FOUND VIA THE LEAST SQUARES
C*                 ANALYSIS.  THEY SHOULD BE DIMENSIONED AT LEAST MC.
C*  AMPC(I),PHG(I) ARE ARRAYS CONTAINING THE AMPLITUDE AND PHASE FOR
C*                 CONSTITUENT NAME(I) AFTER CORRECTIONS FOR NODAL
C*                 MODULATION, ASTRONOMICAL ARGUMENT AND INFERRED
C*                 CONSTITUENTS.  THEIR MINIMUM DIMENSION SHOULD BE MC.
C*  ITH(I),ITM(I), ARE ARRAYS CONTAINING THE TIMES IN HOURS AND MINUTES
C*  ht(I)          AND HEIGHTS, OF THE OBSERVED DATA AS IT IS INPUT BY
C*                 RECORD.  THEY SHOULD BE DIMENSIONED ACCORDINGLY( AT
C*                 PRESENT ONLY 6 OBSERVATIONS ARE EXPECTED PER RECORD).
C*  X(I),Y(I)      ARE ARRAYS CONTAINING ALL THE TIMES(IN HOURS AS
C*                 MEASURED FROM THE CENTRE OF THE ANALYSIS PERIOD) AND
C*                 HEIGHTS OF THE OBSERVED DATA.  THEIR MINIMUM
C*                 DIMENSION SHOULD BE NOBS.
C*  NSTN(I)        IS THE ARRAY CONTAINING THE TIDAL STATION NAME.  IT
C*                 SHOULD HAVE MINIMUM DIMENSION 5.
C*  Q(I)           IS THE OVERDETERMINED ARRAY OF EQUATIONS THAT IS
C*                 SOLVED IN THE LEAST SQUARES SENSE BY THE MODIFIED
C*                 GRAM-SCHMIDT ALGORITHM.  IT SHOULD HAVE THE EXACT
C*                 DIMENSION OF NMAXPM BY NMAXP1.
C*  P(I)           IS THE ARRAY CONTAINING THE TIDAL CONSTITUENT SINE
C*                 AND COSINE COEFICIENTS AS FOUND WITH THE LEAST
C*                 SQUARES FIT.  IT SHOULD HAVE MINIMUM DIMENSION NVAR.
C
C   COMMENT: TO ALIGN CHARACTER*5 ARRAY NAME(MC) IN A COMMON BLOCK IT IS
C            ADVISABLE FOR MC BE MULTIPLE OF 4
C***********************************************************************
C     

      real*8 :: sig_ls(nmaxp1) = 1.d0 ! todo: understandt this sig


      
      real*8 :: residuals(nobs,ndef,nloc)
      real*8, allocatable :: Q(:,:) 

      real*8, allocatable :: cov_by_var(:,:)

      real*8 dif(nobs)
      real*8 sdevtemp
      integer kdtp(nobs),kd(nobs)
	

	integer imin
      
      character*3,dimension(2) :: complabel = (/ '[1]', '[2]' /)
      
      real*8, allocatable :: P(:,:,:)
      real*8 :: cenhr,cumhr,yy

c
c     Additional arrays, for use in the SVD routine (J.Ch., Aug. 1997)
      real*8, allocatable :: matu(:,:),matv(:,:)


      real*8 :: cor(nmaxp1,nmaxp1)
c

      integer ikount,jkount,kkount
      integer idef
      integer iloc
      integer itest
      
      integer :: id1,im1,iy1,id2,im2,iy2,ic1,ic2
      
      integer lm

      integer kd1,kh1,kd2,kh2,khm
      integer istn,isec
      real*8 xlat,xlon

      integer ibin,I
      real*8 fac  
      integer lp,kr1

      real*8 vx, ux, fx

      integer lats, lons,l,k
      real*8 xmid
      integer irep
      integer j,jj1,icode,kh
      integer inflag,kinf
      real*8 arg,arg1,arg2,arg3
      real*8 vxi,uxi
      real*8 fxi,c2,s2,c3,s3
      integer :: nvar         ! number of variables/parameters in overdetermined system
      integer :: meq          ! number of equations in system, generally meq = nobs
      integer :: jcode
      real*8 aamp

      
      integer imax,i2,i21

      integer ii1,i1
      integer jmax,iconst,jconst,ic,iy,im,id,ih
      real*8 var,add
      integer niter,iter,j0
      real*8 cormax,ac,cx,sx,cy,sy,cxpsy,cymsx,apl,amn,epl
      real*8 emin,gpl,gmn,cxmsy,cypsx
      real*8 gpli,gmni
      integer :: mc_used ! number of constituents including z0 and inferred, not satellites

      
      data lp,kr1/6,7/ ! revised

C***********************************************************************
304	format(A80)
      
      istat = 0

      iy1 = itime(1,1)
      im1 = itime(2,1)
      id1 = itime(3,1)
      ic1 = itime(7,1)
      iy2 = itime(1,nobs)
      im2 = itime(2,nobs)
      id2 = itime(3,nobs)
      ic2 = itime(7,nobs)
      do i = 1,nobs
        call GDAY(itime(3,i),itime(2,i),itime(1,i),itime(7,i),kd(i))
      end do


! todo: moved this
c	mf= number of consituents, excluding linear trend. The constant
c	term, Z0 should be first in the list.
c	itrend= 1 if include linear trend
c	itrend= otherwise, no trend
c	ndef=1 if only 1D field to be analysed (eg., elevations)
c	ndef=2 if 2D field: velocity components, EW followed by NS
c	number of unknowns, M, depends on whether we have a linear trend
	if(itrend.eq.1) then
      nvar=2*MF
	else
	nvar=2*MF-1
	end if
      
	ibin=0

      open(unit=lp,file=file_analysis_out,status='unknown',form='formatted')
      !open(unit=11,file=file_svd,status='unknown',form='formatted')
	!open(unit=25,file=file_fitted,status='unknown',form='formatted')
C

      FAC=twopi/360.
      CALL GDAY(ID1,IM1,IY1,IC1,kd1)
      KH1=24*kd1
      CALL GDAY(ID2,IM2,IY2,IC2,kd2)
      KH2=24*(kd2+1)
      KHM=(KH1+KH2)/2
	hrm=khm
      CENHR=DFLOAT((KH2-KH1)/2)
      
      
	xlat=latd+latm/60.
	istn=9999      
c

c ----------------------------------------------------------------------
c
c	read in the astronomical argument information
c 	and re-calculate the constituent frequencies (latitude not used yet)
c
      ! todo: moving here ok? call vt_read_constituent_db(name_constituent_db) 
      DO I=1,MF
        FREQ(I)=FREQC(I)*twopi
      ENDDO
C
C***********************************************************************
C*  DETERMINE THE CENTRAL HOUR OF THE ANALYSIS PERIOD AND SET UP THE
C*  DEPENDENT AND INDEPENDENT VARIABLES, Y AND X.
C
      lats=int(60.*(60.*(xlat-LATD)-LATM))
      IF(ABS(xlat).LT.5.) xlat=SIGN(5.,xlat)
c
c actually, CUMHR=24.d0*(KD-KHM) (check), but keep same notation as before
c
	do i=1,nobs
        cumhr=-cenhr+24.D0*(kd(i)-kd1)
        time(i)=CUMHR+DFLOAT(ITime(4,i))+(DFLOAT(ITime(5,i))
     1             +DFLOAT(ISEC)/60.d0)/60.d0
        x(i)=time(i)-time(1)
        kdtp(i)=kd(i)
      end do

c
C***********************************************************************
C*  SETTING UP THE OVERDETERMINED MATRIX AND SOLVING WITH MODIFIED SVD
C
      if (.not. allocated(Q))then
          allocate(Q(nobs,nvar))
      end if 
      if (.not. allocated(cov_by_var))then
          allocate(cov_by_var(nvar,nvar))
      end if      
      if (.not. allocated(matv))then
          allocate(matu(nobs,nvar))
          allocate(matv(nvar,nvar))
      end if
      if (.not. allocated(p))then
          allocate(p(nvar,ndef,nloc))
      end if 

      irep=0
      xmid=0.5*(x(1)+x(nobs))      
      


      
      Q=0.d0 ! Array set to zero
      DO I=1,nobs
c	if itrend=1, then
c	first 2 parameters are constant and linear trend (per 365 days)
c	fitted as const+trend(t-tmid) where tmid (=xmid) is the middle time
c	of the analysis period (This makes the constant consistent with z0
c	in the old analysis program)
c	If itrend=0 then the second parameter is is associated with the next 
c	constituent
        Q(I,1)=1.
	  if(itrend.eq.1) then
          Q(I,2)=(x(i)-xmid)/(24.*365.)
	  end if
		icode=1
	  kh=24*kdtp(i)+itime(4,i)
		hr=dfloat(kh)+itime(5,i)/60.d0+itime(6,i)/3600.d0
        CALL SETVUF(hr,xlat,mf,name)
        DO J=2,MF
            CALL VUF(name(j),vx,ux,fx) 
c	check to see if this constituent is to be used for inference
			inflag=0
			kinf=0
			if(ninfer_cnstnts.eq.0) go to 500
			do k=1,ninfer_cnstnts
				if(name(j).eq. infer_main_const(k)) then
				inflag=1
				kinf=k
				go to 500
				end if
			end do
500			continue
		if(inflag.eq.0) then
c			ARG=X(I)*FREQ(J)+ux*twopi
			ARG=(vx+ux)*twopi
			if(itrend.eq.1) then
			jkount=2*(J-1)+1
			else
			jkount=2*(J-1)
			end if
			JJ1=jkount+1
			Q(I,jkount)=COS(ARG)*fx
			Q(I,jj1)=SIN(ARG)*fx
		else
			if(itrend.eq.1) then
			jkount=2*(J-1)+1
			else
			jkount=2*(J-1)
			end if
			JJ1=jkount+1
			ARG1=(vx+ux)*twopi
			Q(I,jkount)=COS(ARG1)*fx
			Q(I,JJ1)=SIN(ARG1)*fx
			do 501 kkount=1,ninfer(kinf)
			CALL VUF(infer_names(kinf,kkount),vxi,uxi,fxi) 
c	freq is radians/hr but infer_freq is cycles/hr 
c			ARG1=X(I)*FREQ(J)+ux*twopi
c			ARG2=(X(I)*infer_freq(kinf)+uxi)*twopi
			ARG2=(vxi+uxi)*twopi
			c2=cos(arg2)
			s2=sin(arg2)
			arg3=infer_zeta(kinf,kkount)*fac
			c3=cos(arg3)
			s3=sin(arg3)
			Q(I,jkount)=q(i,jkount)+fxi*infer_amp(kinf,kkount)*(c2*c3-s2*s3)
			Q(I,JJ1)=q(i,jj1)+fxi*infer_amp(kinf,kkount)*(c2*s3+s2*c3)
501			continue
		end if
        ENDDO !j
!		else if(idef.eq.2) then
!		Q(i,nvarp1)=y2(i)
!		icode=2
!		end if ! idef if
      ENDDO  !i

	call log_message( 'assembled overdetermined matrix and/or rhs')
	do iloc=1,nloc
	  do idef=1,ndef
	    sumtest=0.d0
	    do itest=1,nobs
	      sumtest = sumtest+obsdata(itest,idef,iloc)
	      end do
	    end do 
	end do
	
      MEQ=nobs
      av=sum(obsdata,1)/dble(meq)      ! mean of each column of observations
      do iloc=1,nloc
	do idef=1,ndef
        ssq(idef,iloc)=1.0
        resmax_ls(idef,iloc)=1.0
c
C***********************************************************************
C*  CALCULATION OF THE STANDARD DEVIATION OF THE RIGHT HAND SIDES OF
C*  THE OVERDETERMINED SYSTEM
C
        sdevtemp=0.D0
        do I=1,MEQ
          sdevtemp=sdevtemp+(obsdata(i,idef,iloc)-av(idef,iloc))**2.d0
        end do
        sdevtemp=sqrt(sdevtemp/dble(meq-1))
        sdev0(idef,iloc)=sdevtemp
      end do
      end do

C***********************************************************************
C   USE SINGULAR-VALUE-DECOMPOSITION TO SOLVE THE OVERDETERMINED SYSTEM
C
c	do we need to adjust the value for toler ?
c     TOLER=1.d-4
      TOLER=1.d-5
      ! todo: understand this? 
      do i=1,nvar
         sig_ls(i)=1.d0
      enddo


c
c	no solution if meq lt m. ie underdetermined system
c	go to next time series
	if(meq.le.nvar) then
	    call log_message(' underdetermined system: no svd solution')
	    istat = -1
	end if

	call log_message(' applying svd')
      call SVD(Q,matu,matv,cov,eigval,p,obsdata,ndef*nloc,sig_ls,icode,meq,nvar,meq,nvar,
     1        toler,jcode,ssq,resmax_ls,residuals)
      IF(JCODE.GT.0) WRITE(LP,55)JCODE
   55 FORMAT('COLUMN',I5,' IS THE 1ST DEPENDENT COLUMNS IN SVD')
   

c	write out eigenvalues
	wmax=-1000.d0
	wmin=1000.d0
	do i=1,nvar
	  if(eigval(i).gt.wmax) wmax=eigval(i)
	  if(eigval(i).lt.wmin) wmin=eigval(i)
	end do

C***********************************************************************
      do 27 iloc=1,nloc
      do 27 idef=1,ndef
	if(ssq(idef,iloc).gt.1.d-10) then
        rmsr0(idef,iloc)=sqrt(ssq(idef,iloc)/(meq-nvar))
	else
	  rmsr0(idef,iloc)=0.d0
	end if

c	
c***********************************************************************
c	re-calculate the residual sum of squares again just to check svd routine
c	residuals are stored in (q(i,nn),i=1,m)
c

	rmsr(idef,iloc)=0.d0
	resmax(idef,iloc)=0.d0
	iresmax(idef,iloc)=0
      do i=1,nobs
		yy=residuals(i,idef,iloc)
		rmsr(idef,iloc)=rmsr(idef,iloc)+yy*yy
	  if(abs(yy).gt.resmax(idef,iloc)) then
	    resmax(idef,iloc)=abs(yy)
	    imax=i
	  end if
      end do
      iresmax(idef,iloc) = imax
	if(rmsr(idef,iloc).gt.1.d-10) then
	  rmsr(idef,iloc)=sqrt(rmsr(idef,iloc)/dble(nobs-nvar))
	else
	  rmsr(idef,iloc)=0.d0
	end if
c
C***********************************************************************
C*  CALCULATE AMPLITUDES AND PHASES
C
c	if itrend=1 then the linear trend is shown as the phase of the constant 
c	Z0 term (& the true phase of Z0 is zero)
c	otherwise, the phase of Z0 is shown as zero
      AMP(1)=P(1,idef,iloc)
	if(itrend.eq.1) then
      PHI(1)=P(2,idef,iloc)
	else
	PHI(1)=0.
	end if
      do 39 I=2,MF
	  if(itrend.eq.1) then
          I2=2*(I-1)+1
	  else
	    I2=2*(I-1)
	  end if
        I21=I2+1
        C=P(I2,idef,iloc)
        S=P(I21,idef,iloc)
        AAMP=SQRT(C*C+S*S)
        IF(AAMP.LT.1.d-5) GOTO 40
        PHI(I)=ATAN2(S,C)/FAC
        IF(PHI(I).LT.0.) PHI(I)=PHI(I)+360.
        GOTO 39
   40   PHI(I)=0.
   39   AMP(I)=AAMP
C***********************************************************************
c	Note that with f & u included in the lsq fit, we only need V from routine VUF
C	but we don't want to correct with V for a central hour. Better to include
c	the right V in the lsq fit. This has been done.
c      CALL SETVUF(hrm,KON,VX,ux,FX,xlat)
        ampc(1,idef,iloc)=amp(1)
        phg(1,idef,iloc)=phi(1)
      do 45 i=2,mf
c       CALL VUF(hrm,NAME(I),VX,ux,FX,xlat)
        ampc(i,idef,iloc)=amp(i)
	  phg(i,idef,iloc)=phi(i)
c       PHG(I)=VX*360.+PHI(I)
c       PHG(I)=AMOD(PHG(I),360.)
   45 continue
c	write out results for constant term & linear trend
	i=1
	if(cov(1,1).gt.1.d-8) then
	  sig(1,i,idef,iloc)=sqrt(cov(1,1))*rmsr0(idef,iloc)
	else
	  sig(1,i,idef,iloc)=0.d0
	end if
	if(itrend.eq.1.and.cov(2,2).gt.1.d-8) then
	  sig(2,i,idef,iloc)=sqrt(cov(2,2))*rmsr0(idef,iloc)
	else
	  sig(2,i,idef,iloc)=0.d0
	end if
	sig(3,i,idef,iloc)=0.d0
	sig(4,i,idef,iloc)=0.d0
c	results for the other constituents
      do i=2,mf
	  if(itrend.eq.1) then
		  ikount=2*(I-1)+1
	  else
		  ikount=2*(I-1)
	  end if
	  ii1=ikount+1
c	multiply cov values with residual standard deviation, as described in equation
c	(6) of Cherniasky et al. (2001)
        if(cov(ikount,ikount).gt.1.d-8) then
	    sig(1,i,idef,iloc)=sqrt(cov(ikount,ikount))*rmsr0(idef,iloc)
	  else
	    sig(1,i,idef,iloc)=0.d0
	  end if
	  if(cov(ii1,ii1).gt.1.d-8) then
	    sig(2,i,idef,iloc)=sqrt(cov(ii1,ii1))*rmsr0(idef,iloc)
	  else
	    sig(2,i,idef,iloc)=0.
	  end if
c	from equation 11 in Pawlowicz et al (2002)
	  c=ampc(i,idef,iloc)*cos(phg(i,idef,iloc)*fac)
	  s=ampc(i,idef,iloc)*sin(phg(i,idef,iloc)*fac)
        sig(3,i,idef,iloc)=sqrt(((c*sig(1,i,idef,iloc))**2.d0+(s*sig(2,i,idef,iloc))**2.d0)/(c**2+s**2.d0))
        sig(4,i,idef,iloc)=ampc(i,idef,iloc)/sig(3,i,idef,iloc)
      end do
C***********************************************************************
C*    Results for inferred parameters
C
	if(ninfer_cnstnts /= 0)then 
	  l=0
	  do k=1,ninfer_cnstnts
	    do i=2,mf
	      if(name(i).eq.infer_main_const(k)) exit
          end do
          i1=i
	    do kkount=1,ninfer(k)
	      l=l+1
	      ampci(k,kkount,idef,iloc)=ampc(i1,idef,iloc)*infer_amp(k,kkount)
	      phgi(k,kkount,idef,iloc)=phg(i1,idef,iloc)-infer_zeta(k,kkount)
 	    end do 
        end do
      end if
c
C***********************************************************************
C*  RECALCULATE THE RESIDUAL ROOT MEAN SQUARE ERROR
c	using re-constructed time series
C
      ssq(idef,iloc)=0.d0
      
      DO 68 I=1,nobs
        kh=24*kdtp(i)+itime(4,i)
	  hr=24.d0*kdtp(i)+itime(4,i)+itime(5,i)/60.d0+itime(6,i)/3600.d0
        CALL SETVUF(hr,xlat,mf,name)
	if(itrend.eq.1) then
        sum1=p(1,idef,iloc)+p(2,idef,iloc)*(x(i)-xmid)/(365.d0*24.d0)
	else
        sum1=p(1,idef,iloc)
	end if
      do J=2,MF
        CALL VUF(name(j),vx,ux,fx)
c       arg=x(i)*freq(j)+ux*twopi-phi(j)*fac
        arg=(vx+ux)*twopi-phg(j,idef,iloc)*fac
        add=fx*ampc(j,idef,iloc)*cos(arg)
        sum1=sum1+add
      end do
	if(ninfer_cnstnts.eq.0) exit
	do k=1,ninfer_cnstnts
		do kkount=1,ninfer(k)
          CALL VUF(infer_names(k,kkount),vx,ux,fx)
c         arg=(x(i)*infer_freq(k)+ux)*twopi-phgi(k)*fac
          arg=(vx+ux)*twopi-phgi(k,kkount,idef,iloc)*fac
          add=fx*ampci(k,kkount,idef,iloc)*cos(arg)
		  sum1=sum1+add
        end do
      end do
      ic=itime(7,i)
	iy=itime(1,i)
	im=itime(2,i)
	id=itime(3,i)
	ih=itime(4,i)
	imin=itime(5,i)
	dif(i)=obsdata(i,idef,iloc)-sum1
      ssq(idef,iloc)=ssq(idef,iloc)+dif(i)**2
      fitted(i,idef,iloc)=sum1
   68 continue
145	format(6i2,4f10.4)   
146	format(i10,4f10.4)
      rmsr(idef,iloc)=sqrt(ssq(idef,iloc)/(nobs-nvar))
      sdev(idef,iloc)=sqrt(ssq(idef,iloc)/nobs)
      var=ssq(idef,iloc)/nobs
      do j=1,nvar
        do i=1,nvar
          cov_by_var(i,j)=cov_by_var(i,j)*var
        enddo
      enddo

      j0=j0+1
27	continue  ! idef and iloc loops
c
c	compute ellipse parameters if ndef=2
	if(ndef > 1) then
	do iloc=1,nloc
	  do idef=1,ndef
	  if(ampc(1,idef,iloc)>= 0.d0) then
	    phg(1,idef,iloc)=0.d0
	  else
	    phg(1,idef,iloc)=180.d0
	    ampc(1,idef,iloc)=-ampc(1,idef,iloc)
	  end if
	end do !idef

	DO 149 I=1,MF
	cx=ampc(i,1,iloc)*cos(phg(i,1,iloc)*fac)
	sx=ampc(i,1,iloc)*sin(phg(i,1,iloc)*fac)
	cy=ampc(i,2,iloc)*cos(phg(i,2,iloc)*fac)
	sy=ampc(i,2,iloc)*sin(phg(i,2,iloc)*fac)
	cxpsy=0.5d0*(cx+sy)
	cymsx=0.5d0*(cy-sx)
	cxmsy=0.5d0*(cx-sy)
	cypsx=0.5d0*(cy+sx)
	apl=sqrt(cxpsy**2.d0+cymsx**2.d0)
	amn=sqrt(cxmsy**2.d0+cypsx**2.d0)
	if(apl > 1.d-5) then                !todo: were any of these inequalities supposed to be >= or <=?
	  epl=atan2(cymsx,cxpsy)/fac
	else
	  epl=0.d0
	end if
	if(amn > 1.d-5) then
	  emin=atan2(cypsx,cxmsy)/fac
	else
	  emin=0.d0
	end if
	amaj(i)=apl+amn
	amin(i)=apl-amn
	gpl=-epl
	gmn=emin
	if(gmn-gpl.lt.0.) gmn=gmn+360.d0
	g(i)=0.5d0*(gpl+gmn)
	ainc(i)=0.5d0*(gmn-gpl)
	if(gmn < 0.d0) gmn=gmn+360.d0
	if(gmn > 360.d0) gmn=gmn-360.d0
	if(gpl < 0.d0) gpl=gpl+360.d0
	if(g(i) < 0.d0) g(i)=g(i)+360.d0
149	continue
c
	if(ninfer_cnstnts > 0) then
	do 153 k=1,ninfer_cnstnts
	do 154 kkount=1,ninfer(k)
	cx=ampci(k,kkount,1,iloc)*cos(phgi(k,kkount,1,iloc)*fac)
	sx=ampci(k,kkount,1,iloc)*sin(phgi(k,kkount,1,iloc)*fac)
	cy=ampci(k,kkount,2,iloc)*cos(phgi(k,kkount,2,iloc)*fac)
	sy=ampci(k,kkount,2,iloc)*sin(phgi(k,kkount,2,iloc)*fac)
	cxpsy=0.5*(cx+sy)
	cymsx=0.5*(cy-sx)
	cxmsy=0.5*(cx-sy)
	cypsx=0.5*(cy+sx)
	apl=sqrt(cxpsy**2+cymsx**2)
	amn=sqrt(cxmsy**2+cypsx**2)
	if(apl.gt.1.d-5) then
	epl=atan2(cymsx,cxpsy)/fac
	else
	epl=0.
	end if
	if(amn.gt.1.d-5) then
	emin=atan2(cypsx,cxmsy)/fac
	else
	emin=0.
	end if
	amaji(k,kkount)=apl+amn
	amini(k,kkount)=apl-amn
	gpli=-epl
	gmni=emin
	if(gmni-gpli.lt.0.) gmni=gmni+360.
	gi(k,kkount)=0.5*(gpli+gmni)
	ainci(k,kkount)=0.5*(gmni-gpli)
	if(gmni.lt.0.) gmni=gmni+360.
	if(gmni.gt.360.) gmni=gmni-360.
	if(gpli.lt.0.) gpli=gpli+360.
	if(gi(k,kkount).lt.0.) gi(k,kkount)=gi(k,kkount)+360.
154	continue
153	continue
      end if      !  ninfer_cnstnts > 0
	
	end do !iloc
	end if !idef > 0


      mc_used = mf
      do k=1,ninfer_cnstnts
        mc_used=mc_used+ninfer(k)
      end do

      close(lp)
      close(kr1)
      close(10)
      close(13)
      if (allocated(Q)) deallocate(Q)
      if (allocated(p)) deallocate(p)
      if (allocated(cov_by_var)) deallocate(cov_by_var)      
      if (allocated(matv))then
          deallocate(matu)
          deallocate(matv)
      end if
	return
      END
