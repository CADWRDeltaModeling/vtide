
      subroutine write_fit_time_series(filename,ipredtime,fitted,ntime,ndef,nloc)
      implicit none
      integer :: ntime, ndef, nloc
      integer i
      integer :: id,im,iy,ic,ih,imin
      integer iloc, idef
      character*80 filename
      integer ipredtime(7,ntime)
      real*8 fitted(ntime,ndef,nloc)
      open(unit=25,file=filename,status='unknown',form='formatted')

      do i=1,ntime
        iy=ipredtime(1,i)
        im=ipredtime(2,i)
        id=ipredtime(3,i)
        ih=ipredtime(4,i)
        imin=ipredtime(5,i)
        ic=ipredtime(7,i)
        write(25,145) ic,iy,"-",im,"-",id," ",ih,":",imin,
     1       ((fitted(i,idef,iloc),idef=1,ndef),iloc=1,nloc)
      end do
145	format(i2.2,i2.2,1a,i2.2,1a,i2.2,1a,i2.2,1a,i2.2,1000f10.4)
      close(25)
      return
      end subroutine
      
      !> Wrote out a file with simple harmonic estimates for each component and station.
      !> This makes an easy input file for prediction and also is good for comparing stations
      !> since their harmonic estimates end up side-by-side
      ! todo: this doesn't handle itrend
      subroutine vt_write_harmonic_output(file_harmonics,
     &                                    station_id,
     &                                    station_name,
     &                                    latd,
     &                                    latm,
     &                                    lond,
     &                                    lonm,
     &                                    time_zone,
     &                                    nobs, itime,
     &                                    mf,ndef,nloc,itrend,
     &                                    name,freqc,ampc,phg,
     &                                    ninfer_cnstnts,ninfer,infer_names,ampci,phgi
     &                                    )
      implicit none
      integer, parameter :: lm = 1010
      integer, parameter :: max_infer = 80   !todo: hardwire
      integer, parameter :: mcc = 70    ! max allowed consituent estimates      
      
      character*(*), intent(in) :: file_harmonics   !< file to which harmonic fit will be written
      integer, intent(in)       :: station_id       !< numerical id for station
      character*20, intent(in)  :: station_name     !< name of station (max 20 characters)
      integer, intent(in) :: latd                   !< representative latitude, degrees
      integer, intent(in) :: latm                   !< representative latitude, minutes part
      integer, intent(in) :: lond                   !< represenative longitude, degrees (note longitude is not part of the analysis)
      integer, intent(in) :: lonm                   !< representative longitude, minutes part
      character*4,intent(in) :: time_zone     !< time zone. this is not used to adjust phase, just a label
      integer,intent(in) :: nobs             !< number of observations
      integer,intent(in) :: itime(7,nobs)    !< obs times ... year,mon,day,hour,min,sec,century
      integer :: mf            !< number of main frequencies in analysis (excudes inferred constituents and satellites
      integer :: ndef          !< number of components, 1 for a scalar and 2 for a vector like velocity
      integer :: nloc          !< number of individual locations analyzed (stations)
      integer :: itrend
      character*5, intent(in), dimension(*) :: name !< main constituent names 
      real*8, intent(in), dimension(*) :: freqc     !< nominal frequency in cycles/hr (just echoes the input, not really used)  
      integer, intent(in) :: ninfer_cnstnts      ! number of main constituents with inferences
      integer, intent(in) :: ninfer(max_infer)   ! number of inferred constituents

      character*5, intent(in), dimension(max_infer,16) :: infer_names  !< names of inferred constituent for each inference
      real*8, dimension(max_infer,16) :: infer_freq                    !< frequency (c/hr) names of inferred constituent for each inference
      real*8, intent(in) :: ampci(10,10,ndef,nloc)                     !< amplitude of inferred constituents !todo: 10,10 hardwire
      real*8, intent(in) :: phgi(10,10,ndef,nloc)                      !< phase of inferred constituents
      real*8, intent(in) :: ampc(MCC,ndef,nloc)                        !< amplitude of main constituents
      real*8, intent(in) :: phg(MCC,ndef,nloc)                         !< phase of main constituents

      integer  :: iy1,im1,id1,ih1,ic1,iy2,im2,id2,ic2
      

      integer i,jj,kk,k,iloc,idef,kkount
      integer :: mc_used     !<  total constituents to be written, including both main and inferred constituents      
      character*3,dimension(2) :: complabel = (/ '[1]', '[2]' /)
      
      iy1 = itime(1,1)
      im1 = itime(2,1)
      id1 = itime(3,1)
      ih1 = itime(4,1)
      ic1 = itime(7,1)
      iy2 = itime(1,nobs)
      im2 = itime(2,nobs)
      id2 = itime(3,nobs)
      ic2 = itime(7,nobs)      
      
      
      ! calculate total number of constituents including both main and inferred
      mc_used = mf
      do k=1,ninfer_cnstnts
        mc_used=mc_used+ninfer(k)
      end do
      
      open(unit=lm,file=file_harmonics,action='write', status='unknown')
      write(lm,158) station_id,station_name,latd,latm,lond,lonm,mc_used,ndef,nloc,itrend
 158  format('STATION #: ',I,/,'STATION NAME: ',A20,/,'LATITUDE:  '
     2       ,2I3,/,'LONGITUDE: ',I4,I3,/,'# CONSTITUENTS:',i4,/,'# COMPONENTS:',i4,/,
     3       '# LOCATIONS:',i4,/,'TREND:',i4,1x,'(if present, trend coef is output as phase of Z0)')
      write(lm,159)"ANALYSIS START: ",ic1,iy1,"-",im1,"-",id1," ",itime(5,1),":",itime(6,1)
      write(lm,159)"ANALYSIS END:   ",ic2,iy2,"-",im2,"-",id2," ",itime(5,nobs),":",itime(6,nobs)
      write(lm,'(a10,a4)')"TIME ZONE: ",time_zone
 159  format(a16,i2.2,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)
      
      if (ndef > 1) then
        write(lm,'(a5,1x,a11,1x,1000(i4,a3,a5))')"Const","Cycle/hr",
     1   ((iloc,complabel(idef),"Ampl", idef=1,ndef),iloc=1,nloc),
     2   ((iloc,complabel(idef),"Phase", idef=1,ndef),iloc=1,nloc)
      else
        write(lm,'(a5,1x,a11,1x,1000(i4,"_",a7))')"Const","Cycle/hr",
     1   ((iloc,"Ampl   ", idef=1,ndef),iloc=1,nloc),
     2   ((iloc,"Phase  ", idef=1,ndef),iloc=1,nloc)      
      end if
      do i = 1,mf
        write(lm,'(a5,1x,1000f12.7)') name(i),freqc(i),((ampc(i,jj,kk),jj=1,ndef),kk=1,nloc), 
     &    ((phg(i,jj,kk),jj=1,ndef),kk=1,nloc)
      end	do 
      do k=1,ninfer_cnstnts
	  do kkount=1,ninfer(k)
	  write(lm,'(a5,1x,400f12.7)') infer_names(k,kkount),infer_freq(k,kkount),
     1              ampci(k,kkount,1:ndef,1:nloc),
	2              phgi(k,kkount,1:ndef,1:nloc)
        end do
      end do
      close(lm)
      end subroutine      
      

       subroutine vt_write_analysis(filename,
     &                      station_id,
     &                      station_name,       
     &                      latd, latm,  
     &                      lond, lonm, 
     &                      time_zone,
     &                      nobs,itime,
     &                      mf,ndef,nloc,itrend,
     &                      name,freqc,ampc,phg, 
     &                      ninfer_cnstnts,ninfer,
     &                      infer_names,infer_freq,ampci,phgi,
     &                      amaj,amin,ainc,g,
     &                      amaji,amini,ainci,gi,
     &                      eigval,sig,cov,sdev,sdev0,ssq,rmsr,rmsr0,resmax,resmax_ls,iresmax)


      use constituent
      implicit none

      integer, parameter :: mcc = 70         ! todo: hardwire, max allowed consituent estimates
      integer, parameter :: max_infer = 80   !todo: hardwire
      integer, parameter :: nmaxp1 = mc*2
      
      character*(*), intent(in) :: filename  !< name of output file for the write
      integer,intent(in) :: station_id       !< integer id of station
      character*20,intent(in) :: station_name !< name of the station
      integer,intent(in) :: latd              !< representative latitude, degree part
      integer,intent(in) :: latm              !< representative latitude, minute part
      integer,intent(in) :: lond              !< representative longitude, degree part (lon is not used in tidal analysis)
      integer,intent(in) :: lonm              !< representative longitude, minute part
      character*4,intent(in) :: time_zone     !< time zone. this is not used to adjust phase, just a label
      integer,intent(in) :: nobs             !< number of observations
      integer,intent(in) :: itime(7,nobs)    !< obs times ... year,mon,day,hour,min,sec,century   
      
      integer,intent(in) :: mf               !< number of main (ie, not inferred or satellite) frequencies in analysis
      integer,intent(in) :: ndef             !< number of components per station, 1 is scalar and 2 is for u,v pairs. other cases not handled 
      integer,intent(in) :: nloc             !< number of locations or stations
      
      integer,intent(in) :: itrend           !< flag indicating a linear trend should be included: 0 = no, 1 = yes     
      
      character*5,intent(in),dimension(*) :: name    !< name of the main constituents in the analysis
      real*8,intent(in), dimension(*)     :: freqc   !< frequencies of the main constituents. Note these are not actually used because the astronomical argument is computed directly.
      real*8, intent(in) :: ampc(mcc,ndef,nloc)      !< amplitude estimate for each main contituent
      real*8, intent(in) :: phg(mcc,ndef,nloc)       !< phase estimate for each main constituent in degrees      
      
      integer,intent(in)                             :: ninfer_cnstnts    !< number of main constituents that have inferences attached to them
      integer,intent(in)                             :: ninfer(max_infer) !< number of inferences for each main constituent with inferences      
      !character*5, dimension(max_infer),intent(in)   :: infer_main_const !< name of main constituent associated with each inference
      character*5,dimension(max_infer,16),intent(in) :: infer_names      !< name of inferred constituent associated with each inference
      real*8, dimension(max_infer,16),intent(in)     :: infer_freq       !< frequency of inferred constituent associated with each inference
      real*8, intent(in) :: ampci(10,10,ndef,nloc)           !< estimated amplitude of inferred constituents !todo: fix hardwired dimensions
      real*8, intent(in) :: phgi(10,10,ndef,nloc)            !< estimated phase of inferred constituents  
	real*8, intent(in) :: amaj(mcc,nloc)                        !< major ellipse parameter, main constituent            
	real*8, intent(in) :: amin(mcc,nloc)                        !< minor ellipse parameter, main constituent
	real*8, intent(in) :: ainc(mcc,nloc)                        !< incl angle, main constituent
	real*8, intent(in) :: g(mcc,nloc)                           !< phase, main constituent 
	real*8, intent(in) :: amaji(10,10,nloc)                !< major ellipse param, inferred constituent
	real*8, intent(in) :: amini(10,10,nloc)                !< minor ellipse param, inferred constituent
	real*8, intent(in) :: ainci(10,10,nloc)                !< incl angle ellipse param, inferred constituent
	real*8, intent(in) :: gi(10,10,nloc)                   !< phase, inferred constituent
	
      real*8 :: eigval(nmaxp1)     
      real*8,intent(in) ::  sig(4,nmaxp1,ndef,nloc)          !< significance statistics for each coefficient estimate, sig1= ... sig2 = ... sig3 = ... sig4=ttest
      real*8,intent(in) :: cov(2*mf-1+itrend,2*mf-1+itrend)  !< covariance matrix of the variables in the least squares problem (z0 and linear trend if used plus 2 per constituent)

      real*8,intent(in)  :: sdev(ndef,nloc)                   !< standard deviation !todo ... of what?
      real*8,intent(in)  :: sdev0(ndef,nloc)                  !< standard deviation of the original right hand side
      real*8,intent(in) :: ssq(ndef,nloc)                    !< sum squared error of each fit
      real*8,intent(in) :: rmsr(ndef,nloc)                   !< root mean square error of each fit
      real*8,intent(in) :: rmsr0(ndef,nloc)                  !< root mean square of ???
      real*8,intent(in) :: resmax(ndef,nloc)                 !< maximume residual from each rhs
      real*8,intent(in) :: resmax_ls(ndef,nloc)              !< maximume residual from each rhs as given by svd routine (should be same as resmax)
      integer,intent(in):: iresmax(ndef,nloc)               !< index of the largest residual of each fit
      
      
      
      integer :: id1,im1,iy1,id2,ih1,imin1,im2,iy2,ic1,ic2,ih2,imin2
      
      integer, parameter :: lm = 1010 ! output file for array-like list of frequencies and coefs
      integer, parameter :: lp = 1016


C***********************************************************************
	

      integer :: kkount,idef,iloc,i
      integer j,k

      integer, parameter :: nlargest = 20
      real*8 :: xlat, xlon

      integer :: nvar         ! number of variables/parameters in overdetermined system
      integer :: meq          ! number of equations in system, generally meq = nobs


      real*8 :: wmax,wmin
      integer imax
      real*8 :: cor(nmaxp1,nmaxp1)  !todo: hardwired size ... should be just passed in
      integer jmax,iconst,jconst
      integer :: iter
      real*8  cormax



C***********************************************************************
304	format(A80)
      
      ic1 = itime(7,1)
      iy1 = itime(1,1)
      im1 = itime(2,1)
      id1 = itime(3,1)
      ih1 = itime(4,1)
      imin1 = itime(5,1)

      ic2 = itime(7,nobs)      
      iy2 = itime(1,nobs)
      im2 = itime(2,nobs)
      id2 = itime(3,nobs)
      ih2 = itime(4,nobs)
      imin2 = itime(5,nobs)

	xlat=latd+latm/60.
	xlon=lond+lonm/60.

! todo: moved this
c	mf= number of consituents, excluding linear trend. The constant
c	term, Z0 should be first in the list.
c	itrend= 1 if include linear trend
c	itrend= otherwise, no trend
c	ndef=1 if only 1D field to be analysed (eg., elevations)
c	ndef=2 if 2D field: velocity components, EW followed by NS
c	number of unknowns, M, depends on whether we have a linear trend
	if(itrend.eq.1) then
      nvar=2*mf
	else
	nvar=2*mf-1
	end if
      
      open(unit=lp,file=filename,status='unknown',form='formatted')

      write(LP,15) ic1,iy1,im1,id1,ih1,imin1,ic2,iy2,im2,id2,ih2,imin2,time_zone    
 15   format(/'THE ANALYSIS PERIOD IS FROM ',I2.2,I2.2'-',I2.2,'-',I2.2," ",I2.2,":",I2.2,
     1' TO ',I2.2,I2.2,'-',I2.2,'-',I2.2," ",I2.2,":",I2.2,'  IN THE TIME ZONE ',A4)
      write(LP,*)'USING SVD TO SOLVE THE OVERDETERMINED SYSTEM'
c      write(lp,150)ID1,IM1,IC1,IY1,ID2,IM2,IC2,IY2
150	format(2i3,2i2,5x,2i3,2i2)

  
c
c
c ----------------------------------------------------------------------
c
c	read in the astronomical argument information
c 	and re-calculate the constituent frequencies (latitude not used yet)
c
	print *, ' number of observations=',nobs
C
C***********************************************************************
C*  DETERMINE THE CENTRAL HOUR OF THE ANALYSIS PERIOD AND SET UP THE
C*  DEPENDENT AND INDEPENDENT VARIABLES, Y AND X.
C
      write(LP,8) station_id,station_name,latd,latm,lond,lonm
    8 format('STATION # ',I,', ',A20,' LATITUDE',2I3,', LONGITUDE',I4,I3)

      write(LP,255)nobs
  255 format('NUMBER OF POINTS IN THE ANALYSIS =',I6)
      write(lp,*) ' nin=',ninfer_cnstnts     
	
      meq=nobs
      
c	write out eigenvalues
	wmax=-1000.d0
	wmin=1000.d0
	do i=1,nvar
	  if(eigval(i).gt.wmax) wmax=eigval(i)
	  if(eigval(i).lt.wmin) wmin=eigval(i)
	end do
	write(6,*) ' max, min eigenvalues =',wmax,wmin
56	format(10e12.5)
C***********************************************************************

      do 27 iloc=1,nloc
      write(LP,51)iloc
   51 format(//,'*************************',//,'LOCATION #: ',I5,4x,'RESULTS')	
      do 27 idef=1,ndef
      write(LP,52) resmax_ls(idef,iloc),ssq(idef,iloc)
   52 format('LARGEST RESIDUAL MAGNITUDE & RESIDUAL SUM OF SQUARES:'
     1        ,2E12.5)
      write(LP,66) sdev0(idef,iloc),rmsr0(idef,iloc)
   66 format(
     1'ST. DEV. OF RIGHT HAND SIDES OF ORIGINAL OVERDETERMINED SYSTEM:'
     2,E12.5/
     3'                       AND THE ROOT MEAN SQUARE RESIDUAL ERROR:'
     4,E12.5)
c	
c***********************************************************************
c	re-calculate the residual sum of squares again just to check svd routine
c	residuals are stored in (q(i,nn),i=1,m)
c

	write(lp,*) ' rms residual: brute force =',rmsr(idef,iloc)
	write(lp,*) ' max residual: ',resmax(idef,iloc),iresmax(idef,iloc)

        write(lp,41)
   41  format('HARMONIC ANALYSIS RESULTS: AMPLITUDES, PHASE LAGS, C, S, 
     1& amp SD estimates, t-test value')

c	results for main constituents
        do 42 I=1,MF
   42   write(lp,43) name(i),freqc(i),ampc(i,idef,iloc),phg(i,idef,iloc),
     1               sig(1,i,idef,iloc),sig(2,i,idef,iloc),sig(3,i,idef,iloc),sig(4,i,idef,iloc)
   43   format(5x,a5,4x,f12.9,2x,f10.5,2x,f10.3,5x,4f8.3)

	if(ninfer_cnstnts /= 0) then
	  write(lp,*) ' INFERENCE RESULTS'
	    do 75 k=1,ninfer_cnstnts
	      do 751 kkount=1,ninfer(k)
	        write(lp,79) infer_names(k,kkount),infer_freq(k,kkount),ampci(k,kkount,idef,iloc),
	1                     phgi(k,kkount,idef,iloc)
  79	        format(5x,a5,4x,f12.9,15x,f10.4,5x,f10.4)
  751	      continue
  75	    continue
      end if
c
c***********************************************************************
c	compute (Cherniawsky et al (2001), page 653) and rank correlation coefficients
c	largest niter value are computed and shown
c	if itrend=1, then the second part of Z0 is the linear trend coefficient
c
      imax=0
      jmax=0
	do i=1,nvar
	  do j=1,i
	    cor(i,j)=cov(i,j)/sqrt(cov(i,i)*cov(j,j))
	  end do
	end do
	do iter=1,nlargest
	  cormax=0.d0
	  do i=2,nvar
	  do j=1,(i-1)
	    if(abs(cor(i,j))> cormax) then
	      cormax=abs(cor(i,j))
	      imax=i
	      jmax=j
	    end if
	  end do
        end do
	  iconst=(imax+2-itrend)/2
	  jconst=(jmax+2-itrend)/2
	  write(lp,83) iter,cormax,imax,jmax,name(iconst),name(jconst)
83	  format(i5,' largest correlation coefficient is ',f8.3,' at (i,j)='
     1,2i5,' for constituents ',a5,' and ',a5) 
	  cor(imax,jmax)=0.d0  !todo: this destroys cor   
      end do

c
C***********************************************************************
C*  RECALCULATE THE RESIDUAL ROOT MEAN SQUARE ERROR
c	using re-constructed time series
C
      write(LP,70)nobs,nvar,'    ',xlat,xlon,sdev0(idef,iloc),sdev(idef,iloc)
   70   format('N,M,LAT,LON,SDEV0,SDEV:  ',2i10,a4,f9.4,f10.4,2f10.2)
        if(ninfer_cnstnts>0) then
        write(LP,71) rmsr(idef,iloc)
   71   format('ROOT MEAN SQUARE RESIDUAL ERROR AFTER INFERENCE IS',
     1          E15.6, //)
        else
        write(LP,72) rmsr(idef,iloc)
   72   format('RECALCULATED ROOT MEAN SQUARE RESIDUAL ERROR IS   ',
     1          E15.6, //)
        endif
c
27	continue  ! idef and iloc loops
c
c	compute ellipse parameters if ndef=2
	if(ndef > 1) then
	do iloc=1,nloc
	write(lp,*) ' ELLIPSE PARAMETERS:                        Major, 
	1  Minor, angle incl,  phase'
	write(lp,*) ' Constituents included directly (not inferred)' 
	do 149 I=1,MF
      write(LP,151) name(i),freqc(i),amaj(i,iloc),amin(i,iloc),ainc(i,iloc),g(i,iloc)
151   format(5X,A5,4X,F12.9,15X,2F10.4,2F10.2)
149	continue
c
	if(ninfer_cnstnts > 0) then
	write(lp,*) ' Inferred constituents' 
	do 153 k=1,ninfer_cnstnts
	do 154 kkount=1,ninfer(k)
        write(LP,151) infer_names(k,kkount),infer_freq(k,kkount),
     1  amaji(k,kkount,iloc),amini(k,kkount,iloc),
	1  ainci(k,kkount,iloc),gi(k,kkount,iloc)
154	continue
153	continue
      end if      !  ninfer_cnstnts > 0
	
	end do !iloc
	end if !idef > 0

      close(lp)
	return
      end subroutine
