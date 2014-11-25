      module file_io
      character*80 ::  file_analysis_out !  6
      character*80 ::  file_constituent_out !  6      
      character*80 ::  file_analysis     !  7
      character*80 ::  file_const_db     !  8
      character*80 ::  file_obs          !  9
      character*80 ::  file_fitted_out       ! 25
      integer :: verbosity = 0
      
      contains
      
      subroutine log_message(message)
      implicit none
      character*(*) :: message !< message to be logged
      if (verbosity > 0)then
          print*,message
      end if
      return
      end subroutine 



C***********************************************************************
C*  READ FROM DEVICE KR1 THE ANALYSIS TYPE AND TIDAL STATION DETAILS.
C*
C*  1)ONE RECORD FOR THE VARIABLES MF 
C*  MF   = THE NUMBER OF CONSTITUENTS, INCLUDING THE CONSTANT TERM Z0,
C*         TO BE IN THE LEAST SQUARES FIT.
C*
C*  2)ONE RECORD FOR EACH OF THE MF CONSTITUENTS TO BE INCLUDED IN THE
C*  FIT.  EACH RECORD CONTAINS THE VARIABLES NAME AND FREQC IN THE
C*  FORMAT (A5,2X,F13.10).  NAME IS THE CONSTITUENT NAME, WHICH SHOULD
C*  BE LEFT JUSTIFIED IN THE ALPHANUMERIC FIELD, WHILE FREQC IS ITS
C*  FREQUENCY MEASURED IN CYCLES PER HOUR.
C*
C*  3) ONE RECORD IN THE FORMAT (8I5) CONTAINING THE FOLLOWING
C*  INFORMATION ON THE TIME PERIOD OF THE ANALYSIS.
C*  ID1,IM1,IY1 - DAY,MONTH,YEAR OF THE BEGINNING OF THE ANALYSIS
C*                PERIOD,
C*  ID2,IM2,IY2 - DAY,MONTH,YEAR OF THE END OF THE ANALYSIS PERIOD.
C*  IC1,IC2     - CENTURY OFR THE BEGINNING AND END OF THE ANALYSIS
C*                PERIOD (ZERO VALUES ARE RESET TO 19)
C*
C*
C*  4)ONE RECORD IN THE FORMAT (I5,5A4,1X,A4,4I5) CONTAINING THE
C*  FOLLOWING TIDAL STATION  INFORMATION.
C*  JSTN      = TIDAL STATION NUMBER,
C*  (NSTN(I),I=1,5) = TIDAL STATION NAME,
C*  ITZ       = TIME ZONE IN WHICH THE OBSERVATIONS WERE RECORDED,
C*  LATD,LATM = STATION LATITUDE IN DEGREES AND MINUTES,
C*  LOND,LONM = STATION LONGITUDE IN DEGREES AND MINUTES.
C*
C*  5)ONE SET RECORDS FOR EACH POSSIBLE INFERENCE. THE FIRST RECORD HAS THE 
C*	CONSITUENT NAME, ITS FREQUENCY, AND THE NUMBER OF CONSTITUENTS TO BE
C*	INFERRED (4X,A5,E16.10,i5), WHILE THERE IS ONE RECORD FOR EACH OF THE 
C*	CONSTITUENTS TO BE INFERRED WITH THE NAME, FREQUENCY, AMPLITUDE RATIO 
C*	(INFERRED TO REFERENCE) AND PHASE DIFFERENCE (GREENWICH PHASE LAG OF 
C*	THE INFERRED CONSTITUENT SUBTRACTED FROM THE GREENWICH PHASE LAG OF THE
C*    (ANALYSED CONSTITUENT IN THE FORMAT(4X,A5,E16.10,2F10.3)
C*
C*  FOR KR1 INPUT, ALL CONSTITUENT NAMES SHOULD BE LEFT JUSTIFIED IN
C*  THE ALPHANUMERIC FIELD, FREQUENCIES ARE MEASURED IN CYCLES/HOUR, AND
C*  ALL CONSTITUENTS MUST BE INCLUDED IN THE LIST IN READ FROM FNAM8.
C      

      subroutine vt_analysis_info(file_analysis, 
     &                            lat_degree, lat_min,  
     &                            lon_degree, lon_min, 
     &                            station_id, station_name,
     &                            mf,ndef,itrend,nloc,
     &                            cnstnt_names,
     &                            freq_cycles,
     &                            time_zone,
     &                            nin,
     &                            ninf,
     &                            infer_main_const,
     &                            infer_names,
     &                            infer_freq,     
     &                            infer_amp, !todo: was 'R'
     &                            infer_zeta,  !todo: was 'Zeta'
     &                            id1,im1,iy1,id2,im2,iy2,ic1,ic2)
      !use analysis_request
      implicit none
      integer :: nin
      integer, parameter   :: KR1 = 10
      real*8, dimension(*) :: freq_cycles
      real*8, dimension(80,16) :: infer_amp         !todo: hardwired size for these
      real*8, dimension(80,16) :: infer_zeta
      real*8, dimension(80,16) :: infer_freq        ! legacy: used in report, but not in calcs
      real*8, dimension(80) :: infer_main_freq      ! legacy: not used for anything
      integer :: ninf(*)
      character*80 :: file_analysis
      character*5, dimension(*) :: cnstnt_names
      character*5, dimension(*) :: infer_main_const

      character*5, dimension(1000,*) :: infer_names  !todo: dimension correctly
      character*4 :: time_zone
      character*20 :: station_name
      character*1 :: sep1,sep2
      integer :: mf
      integer :: ndef
      integer :: station_id
      integer :: lat_degree
      integer :: lat_min
      integer :: lon_degree
      integer :: lon_min
      integer :: i
      integer :: itrend
      integer :: id1,im1,iy1,id2,im2,iy2,ic1,ic2
      integer :: k,kkount
      integer :: nloc
      open(unit=KR1, file=file_analysis, status='old',action='read')
      read(KR1,*) mf,ndef,itrend,nloc
      print*,"Number of locations is: ", int(nloc)
	! todo: all this reading stuff should be put in a logging function
	!print *, ' number of constituents & degrees of freedom=',mf,ndef
	!if(itrend.eq.1) then
	!print *, ' a linear trend is included in the analysis'
	!else
	!print *, ' no linear trend is included'
	!end if
   10 format(2I5,F5.2)
      read(KR1,*) (cnstnt_names(I),freq_cycles(I),I=1,mf)
      !print *,(cnstnt_names(I),freq_cycles(I),i=1,mf)
   11 format(4x,A5,F16.10)
      read(KR1,7) ic1,iy1,sep1,im1,sep2,id1
      read(KR1,7) ic2,iy2,sep1,im2,sep2,id2
    7 format(i2,i2,1a,i2,1a,i2)
      !print *, id1,im1,iy1,id2,im2,iy2,ic1,ic2
      if(ic1.eq.0) ic1=19
      if(ic2.eq.0) ic2=19

      read(kr1,*) station_id,station_name,time_zone,
     &            lat_degree,lat_min,lon_degree,lon_min
    9 format(I5,1X,A20,1X,A4,4I5)

c	
c	read in inference information now as it will be used in the lsq matrix
c
      DO 1020 K=1,10
      READ(KR1,*,END=1030)infer_main_const(K),infer_main_freq(K),ninf(k)
      IF(len_trim(infer_main_const(K)) == 0) GO TO 1040
	do 1021 kkount=1,ninf(k)
 1021	read(kr1,*) infer_names(K,kkount),infer_freq(k,kkount),
     & infer_amp(K,kkount),infer_zeta(K,kkount)

1020  CONTINUE
1030  continue
1040  NIN=K-1
	!write(lp,*) ' nin=',nin   !todo: restore

      
      return
      end subroutine
      
c**********************************************************************      
      subroutine vt_count_data(file_obs,ndef,nloc,nobs,
     &                         id1,im1,iy1,ic1,
     &                         id2,im2,iy2,ic2)
      implicit none
      character*80, intent(in) :: file_obs
      integer, intent(in) :: ndef
      integer, intent(in) :: id1,im1,iy1,ic1
      integer, intent(in) :: id2,im2,iy2,ic2
      character*16 :: datestr
      character*16000 :: cbuffer      
      character*1 :: sep1,sep2,sep3,sep4      
      integer :: nobs
      integer :: nloc
      integer i,j,k
      integer kd1,kd2
      integer :: kr2 = 32
      integer idd,imm,icc,iyy,ihh,imin,isec,kdd
      real*8 :: observed(ndef,nloc)
      integer :: lp = 6  !todo: hardwired
      open(unit=kr2, file=file_obs)

C***********************************************************************
C*  Count observations

      CALL GDAY(ID1,IM1,IY1,IC1,kd1)
      CALL GDAY(ID2,IM2,IY2,IC2,kd2)
	  i=0
	  kdd = -1
      do while (kdd <= kd2)  !todo: what if we get nowhere close?
      !todo: hardwired max of 400 values
14    read(kr2,'(a16,a)',end=143) datestr,cbuffer
      read(cbuffer,*)((observed(j,k),j=1,ndef),k=1,nloc)
      read(datestr,'(i2,i2,a1,i2,a1,i2,a1,i2,a1,i2)')
     & icc,iyy,sep1,imm,sep2,idd,sep3,ihh,sep4,imin

	  isec=0
      call GDAY(idd,imm,iyy,icc,kdd)
      if(kdd.lt.kd1) then
        cycle ! before analysis period
      end if
      if(kdd.gt.kd2) then
        exit  ! after the analysis period
      end if
	i=i+1
	end do
	
143   nobs=i
	if(nobs.le.10) then
	write(lp,*) ' not enough observations: nobs',nobs
        WRITE(LP,*) icc,iyy,imm,idd,ihh,imin
        WRITE(LP,*)'kd, kd1, kd2 =',kdd,kd1,kd2
	  stop 1
	end if
	close(kr2)
	return
	end subroutine
      
      
      ! todo: be able to directly read data from file that has a relative time in first colunm
      ! by using the start time of the analysis
      subroutine vt_read_data(file_obs,ndef,nloc,
     &                        id1,im1,iy1,ic1,
     &                        id2,im2,iy2,ic2,  
     &                        nobs,itime,obsdata
     &                        )
      implicit none
      character*80, intent(in) :: file_obs
      integer, intent(in) :: ndef
      integer, intent(in) :: nloc
      integer, intent(in) :: id1,im1,iy1,ic1
      integer, intent(in) :: id2,im2,iy2,ic2
      character*16 :: datestr
      character*1 :: sep1,sep2,sep3,sep4
      character*16000 :: cbuffer      
!c      integer, intent(in) :: max_obs
      integer :: itime(7, nobs)
      real*8 obsdata(nobs,ndef,nloc)
      integer :: nobs
      integer i,j,k
      integer kd1,kd2
      integer :: kr2 = 32
      integer idd,imm,icc,iyy,ihh,imin,isec,kdd
      real*8 observed(ndef,nloc)
      integer :: lp = 6  !todo: hardwired
      open(unit=kr2, file=file_obs)

C***********************************************************************
C*  READ FROM DEVICE KR2 THE OBSERVED TIDAL HEIGHTS/velocities AND TIMES.
C*
      CALL GDAY(ID1,IM1,IY1,IC1,kd1)
      CALL GDAY(ID2,IM2,IY2,IC2,kd2)
	i=0
	kdd = -1
      do while (kdd <= kd2)	
14	  read(kr2,'(a16,a)',end=143) datestr,cbuffer
        read(cbuffer,*)((observed(j,k),j=1,ndef),k=1,nloc)
145	  format(a16,1000f10.4)   !todo: hardwired 1000
        read(datestr,'(i2,i2,a1,i2,a1,i2,a1,i2,a1,i2)')
     &   icc,iyy,sep1,imm,sep2,idd,sep3,ihh,sep4,imin
	isec=0      
      call GDAY(idd,imm,iyy,icc,kdd)
      if(kdd.lt.kd1)then
c      if(kdd.lt.kd1.or.kdd.gt.kd2) then
        !WRITE(LP,*) icc,iyy,imm,idd,ihh,imin 
        !WRITE(LP,*)'kd, kd1, kd2 =',kdd,kd1,kd2   ! todo: disabled logging
	!write(lp,*) ' observation before analysis period'
        cycle
      endif
      if (kdd .gt. kd2) exit
	i=i+1
	itime(1,i)=iyy
      itime(2,i)=imm
      itime(3,i)=idd
      itime(4,i)=ihh
      itime(5,i)=imin
      itime(6,i)=isec
	itime(7,i)=icc
	obsdata(i,:,:)=observed
	end do
143	nobs = i
      if(nobs.le.10) then
	write(lp,*) ' not enough observations: nobs',nobs
        WRITE(LP,*) icc,iyy,imm,idd,ihh,imin
        WRITE(LP,*)'kd, kd1, kd2 =',kdd,kd1,kd2
	  stop 1
	end if
	close(kr2)
	return
	end subroutine



	subroutine vt_read_constituent_file_header(fname, 
     1                         station_no, 
     1                         station_name, 
     1                         latd, latm, 
     1                         lond,lonm,
     1                         istart,
     1                         iend,
     1                         nconstituent,ndef,nloc,itrend)
      implicit none
      ! args
      character*80, intent(in) :: fname
      integer, intent(out) :: station_no
      character*20, intent(out) :: station_name
      integer, intent(out) :: latd,latm
      integer, intent(out) :: lond,lonm
      integer, intent(out) :: istart(7)
      integer, intent(out) :: iend(7)
      integer :: nconstituent
      integer, parameter :: max_cnstnt = 80  ! todo: hardwired
      integer, intent(out) :: ndef
      integer, intent(out) :: nloc
      integer, intent(out) :: itrend      
      character*5, dimension(1) :: all_names
      ! locals
      character*128 cbuffer1,cbuffer2
      character*1 :: sep1,sep2,sep3,sep4
      integer, parameter :: lcf = 58
      integer :: i,j,k
      integer :: ic1,iy1,im1,id1,ic2,iy2,im2,id2,ih1,imin1,ih2,imin2
      open(unit = lcf, file=fname, action='read', status='old')
      read(lcf,*)cbuffer1,cbuffer2,station_no
      read(lcf,*)cbuffer1,cbuffer2,station_name
      read(lcf,*)cbuffer1,latd,latm
      read(lcf,*)cbuffer1,lond,lonm
      read(lcf,*)cbuffer1,cbuffer2,nconstituent
      read(lcf,*)cbuffer1,cbuffer2,ndef
      read(lcf,*)cbuffer1,cbuffer2,nloc
      read(lcf,*)cbuffer1,itrend
      read(lcf,722)cbuffer1,ic1,iy1,sep1,im1,sep2,id1,sep3,ih1,sep4,imin1
  722 format(a16,i2,i2,a1,i2,a1,i2,a1,i2,a1,i2)
      istart(1)=iy1
      istart(2)=im1
      istart(3)=id1
      istart(4)=ih1
      istart(5)=imin1
      istart(6)=0
	istart(7)=ic1
      read(lcf,722)cbuffer1,ic2,iy2,sep1,im2,sep2,id2,sep3,ih2,sep4,imin2
	iend(1)=iy2    ! convert hardwired time indices to constants iend(YR_NDX)
      iend(2)=im2
      iend(3)=id2
      iend(4)=ih2
      iend(5)=imin2
      iend(6)=0
	iend(7)=ic2
      read(lcf,*)cbuffer1 ! time zone
      read(lcf,*) cbuffer1
      if (cbuffer1(1:4) /= 'Cons') then
      print*,"Something wrong in header"
      end if
      close(lcf)
      return
      end subroutine
      

	
	subroutine vt_read_constituent_file(fname, 
     1                         station_no, 
     1                         station_name, 
     1                         latd, latm, 
     1                         lond,lonm,
     1                         istart,
     1                         iend,
     1                         nconstituent,ndef,nloc,itrend,
     1                         all_names,
     1                         all_freq,
     1                         all_amp,
     1                         all_phase)
      implicit none
      ! args
      character*80, intent(in) :: fname
      integer, intent(out) :: station_no
      character*20, intent(out) :: station_name
      integer, intent(out) :: latd,latm
      integer, intent(out) :: lond,lonm
      integer, intent(out) :: istart(7)
      integer, intent(out) :: iend(7)
      integer :: nconstituent
      integer, parameter :: max_cnstnt = 80  ! todo: hardwired
      integer, intent(inout) :: ndef
      integer, intent(inout) :: nloc
      integer, intent(inout) :: itrend      
      character*5, dimension(max_cnstnt) :: all_names
      real*8 :: all_freq(max_cnstnt)
      real*8 :: all_amp(max_cnstnt,ndef,nloc)
      real*8 :: all_phase(max_cnstnt,ndef,nloc)
      ! locals
      character*128 cbuffer1,cbuffer2
      character*1 :: sep1,sep2,sep3,sep4
      integer, parameter :: lcf = 58
      integer :: i,j,k
      integer :: ic1,iy1,im1,id1,ic2,iy2,im2,id2,ih1,imin1,ih2,imin2
      
      
      open(unit = lcf, file=fname, action='read', status='old')
      read(lcf,*)cbuffer1,cbuffer2,station_no
      read(lcf,*)cbuffer1,cbuffer2,station_name
      read(lcf,*)cbuffer1,latd,latm
      read(lcf,*)cbuffer1,lond,lonm
      read(lcf,*)cbuffer1,cbuffer2,nconstituent
      read(lcf,*)cbuffer1,cbuffer2,ndef
      read(lcf,*)cbuffer1,cbuffer2,nloc
      read(lcf,*)cbuffer1,itrend
      read(lcf,723)cbuffer1,ic1,iy1,sep1,im1,sep2,id1,sep3,ih1,sep4,imin1
  723 format(a16,i2,i2,a1,i2,a1,i2,a1,i2,a1,i2)
      istart(1)=iy1
      istart(2)=im1
      istart(3)=id1
      istart(4)=ih1
      istart(5)=imin1
      istart(6)=0
	istart(7)=ic1
      read(lcf,723)cbuffer1,ic2,iy2,sep1,im2,sep2,id2,sep3,ih2,sep4,imin2
	iend(1)=iy2    ! convert hardwired time indices to constants iend(YR_NDX)
      iend(2)=im2
      iend(3)=id2
      iend(4)=ih2
      iend(5)=imin2
      iend(6)=0
	iend(7)=ic2
      read(lcf,*)cbuffer1 ! time zone
      read(lcf,*) cbuffer1
      if (cbuffer1(1:4) /= 'Cons') then
      print*,"Something wrong in header"
      end if
      
      do i = 1,nconstituent
      read(lcf,*) all_names(i),all_freq(i),
     1        ((all_amp(i,j,k),j=1,ndef),k=1,nloc),
     1        ((all_phase(i,j,k),j=1,ndef),k=1,nloc)
      end do
      close(lcf)
      return
      end subroutine

   
      end module




