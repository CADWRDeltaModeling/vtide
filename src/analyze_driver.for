
!****************************************************************************
!
!  PROGRAM: versatile_analyze
!
!  PURPOSE:  Main shell program to analyze tides using the versatile system devised
!             by Foreman (***).
!
!****************************************************************************

      program versatile_analyze
      use file_io
      use constituent
	implicit none
	integer , parameter:: kin = 4      !< unit number for the main input file
      character(LEN=80)  :: main_file    
      character*4        :: time_zone
      integer, parameter :: max_main_cnstnt = 80
      integer, parameter :: max_infer = 80
      
      integer, parameter :: mcc = 70       !< max allowed consituent estimates
      integer, parameter :: nmaxp1 = mc*2  !< dimension 'at least as great as the number of variables'
    
      
      
      character*5,dimension(max_infer,16) :: infer_names
      real*8, dimension(max_infer,16)     :: infer_amp
      real*8, dimension(max_infer,16)     :: infer_zeta
      real*8,dimension(max_infer,16)      :: infer_freq
      character*5, dimension(max_infer)   :: infer_main_const  !todo konan, not dimed correctly
      integer :: ninfer_cnstnts
      integer, dimension(max_infer)       :: ninfer
      character*5, dimension(max_main_cnstnt) :: cnstnt_names = ' '
      real*8, dimension(max_main_cnstnt) :: freq_cycles  
      integer :: lat_degree
      integer :: lat_min
      integer :: lon_degree
      integer :: lon_min
      integer :: station_id
      character*20 :: station_name
      integer :: mf
      integer :: ndef
      integer :: nloc
      integer :: itrend
      integer :: nvar ! number of variables in overdetermined system = 2*mf-1+itrend
      integer :: id1,im1,iy1,id2,im2,iy2,ic1,ic2
      integer :: nobs
      integer :: istat
!c      integer, parameter ::max_obs = 106000
      integer, allocatable :: itime(:,:)
      real*8, allocatable :: obsdata(:,:,:)

      real*8,allocatable :: cov(:,:)
      real*8 :: eigval(nmaxp1)
      
      real*8,allocatable :: sig(:,:,:,:)
      real*8,allocatable :: ampci(:,:,:,:)
      real*8,allocatable :: phgi(:,:,:,:)
      real*8,allocatable :: ampc(:,:,:)
      real*8,allocatable :: phg(:,:,:)  
      real*8,allocatable :: sdev(:,:)
      real*8,allocatable :: sdev0(:,:)
      real*8,allocatable :: rmsr(:,:)     
      real*8,allocatable :: rmsr0(:,:)
      real*8,allocatable :: resmax(:,:)     ! max residual recalculated brute force
      real*8,allocatable :: resmax_ls(:,:)  ! max residual as given by least square solver   
      integer,allocatable :: iresmax(:,:)   
      real*8,allocatable :: ssq(:,:)
      real*8,allocatable :: amaj(:,:)
	real*8,allocatable :: amin(:,:)
	real*8,allocatable :: ainc(:,:)
	real*8,allocatable :: g(:,:)
	real*8,allocatable :: amaji(:,:,:)
	real*8,allocatable :: amini(:,:,:)
	real*8,allocatable :: ainci(:,:,:)
	real*8,allocatable :: gi(:,:,:)      
	real*8,allocatable :: fitted(:,:,:) 

      
!c	FILE I/O
!c	KIN is the master input file. 
!c	fnam6 is the file to which the output is sent. It is assigned the number
!c		lp.
!c	fnam7 is file containing the constituents to be included in the analysis,
!c		the analysis period, inference parameters, the flag controlling
!c		height or current analyses, and site information. It is assigned the
!c		number kr1.
!c	fnam8 is the file containing all the astronomical argument information 
!c		(it should not have to be changed)
!c	fnam9 is the file containing the observations and their times. It is
!c		assigned the number kr2. 
!c	fnam11 is a file to which information on the SVD matrix fit is output when
!c		ibin > 0. This info includes matrix covariances which are useful
!c		in determining (in)dependence of the chosen constituents.
!c	Original, fitted, and residual time series are output to file 25.
!c
      call parse_command_line(main_file)
     
      open(UNIT=KIN,FILE=main_file,STATUS='OLD',ACTION='READ')    
      read(KIN,*) file_analysis
      read(KIN,*) file_const_db
      read(KIN,*) file_obs
      read(KIN,*) file_analysis_out
      read(KIN,*) file_constituent_out
      read(KIN,*) file_fitted_out     
      close(KIN)

      call vt_read_constituent_db(file_const_db)

      
      call vt_analysis_info(file_analysis, 
     &                      lat_degree, lat_min,  
     &                      lon_degree, lon_min, 
     &                      station_id, station_name,
     &                      mf,ndef,itrend,nloc,
     &                      cnstnt_names,
     &                      freq_cycles,
     &                      time_zone,
     &                      ninfer_cnstnts,
     &                      ninfer,
     &                      infer_main_const,
     &                      infer_names, 
     &                      infer_freq,      
     &                      infer_amp, 
     &                      infer_zeta,
     &                      id1,im1,iy1,id2,im2,iy2,ic1,ic2)
     
      nvar = 2*mf - 1 + itrend
      print*,"Read analysis (station, constituent, time) specifications"
	if (.not. allocated(cov))   allocate(cov(nvar,nvar))
      if (.not. allocated(sig))    allocate( sig(4,nmaxp1,ndef,nloc))
      if (.not. allocated(ampci))  allocate(ampci(10,10,ndef,nloc))
      if (.not. allocated(phgi))   allocate(phgi(10,10,ndef,nloc))
      if (.not. allocated(ampc))   allocate(ampc(mcc,ndef,nloc))
      if (.not. allocated(phg))    allocate(phg(mcc,ndef,nloc))
      if (.not. allocated(sdev))   allocate(sdev(ndef,nloc))
      if (.not. allocated(sdev0))  allocate(sdev0(ndef,nloc))
      if (.not. allocated(rmsr))   allocate(rmsr(ndef,nloc))
      if (.not. allocated(rmsr0))  allocate(rmsr0(ndef,nloc))
      if (.not. allocated(resmax)) allocate(resmax(ndef,nloc))
      if (.not. allocated(resmax_ls)) allocate(resmax_ls(ndef,nloc))
      if (.not. allocated(iresmax)) allocate(iresmax(ndef,nloc))      
      if (.not. allocated(ssq))    allocate(ssq(ndef,nloc))
      if (.not. allocated(amaj))   allocate(amaj(mcc,nloc))
	if (.not. allocated(amin))   allocate(amin(mcc,nloc))
	if (.not. allocated(ainc))   allocate(ainc(mcc,nloc))
	if (.not. allocated(g))      allocate(g(mcc,nloc))
	if (.not. allocated(amaji))  allocate(amaji(10,10,nloc))
	if (.not. allocated(amini))  allocate(amini(10,10,nloc))
	if (.not. allocated(ainci))  allocate(ainci(10,10,nloc))
	if (.not. allocated(gi)) allocate(gi(10,10,nloc))
      
      call vt_count_data(file_obs,ndef,nloc,nobs,
     &                  id1,im1,iy1,ic1,
     &                  id2,im2,iy2,ic2)
      if (.not. allocated(obsdata)) allocate(obsdata(nobs,ndef,nloc))
      if (.not. allocated(fitted)) allocate(fitted(nobs,ndef,nloc))
      if (.not. allocated(itime)) allocate(itime(7,nobs))      
      call vt_read_data(file_obs,ndef,nloc,
     &                  id1,im1,iy1,ic1,
     &                  id2,im2,iy2,ic2,
     &                  nobs,itime,obsdata)
      
      print*,"Read data complete"
      call vt_analyze(nobs,itime,obsdata,
     &                lat_degree, lat_min,
     &                mf,ndef,itrend,nloc,
     &                cnstnt_names,
     &                freq_cycles,
     &                ninfer_cnstnts,
     &                ninfer,
     &                infer_main_const,
     &                infer_names, 
     &                infer_freq,     
     &                infer_amp, 
     &                infer_zeta,
     &                ampc,phg,
     &                ampci,phgi,
     &                amaj,amin,ainc,g,
     &                amaji,amini,ainci,gi,
     &                eigval,sig,cov,sdev,sdev0,ssq,rmsr,rmsr0,resmax,resmax_ls,iresmax,
     &                fitted,
     &                istat)


      print*,"Analysis complete, writing output"
 
      print*,"Harmonic output"
      call vt_write_harmonic_output(file_constituent_out,
     &                              station_id,
     &                              station_name,
     &                              lat_degree,lat_min,lon_degree,lon_min,time_zone,
     &                              nobs,itime,
     &                              mf,ndef,nloc,itrend,
     &                              cnstnt_names,freq_cycles,ampc,phg,
     &                              ninfer_cnstnts,ninfer,infer_names,
     &                              ampci,phgi) 

       print*,"Analysis output"
       call vt_write_analysis(file_analysis_out, !"out.txt",
     &                         station_id,
     &                         station_name,       
     &                         lat_degree, lat_min,  
     &                         lon_degree, lon_min, 
     &                         time_zone,
     &                         nobs,itime,
     &                         mf,ndef,nloc,itrend,
     &                         cnstnt_names,freq_cycles,ampc,phg, 
     &                         ninfer_cnstnts,ninfer,
     &                         infer_names,infer_freq,
     &                         ampci,phgi,
     &                         amaj,amin,ainc,g,
     &                         amaji,amini,ainci,gi,
     &                         eigval,sig,cov,sdev,sdev0,ssq,rmsr,rmsr0,resmax,resmax_ls,iresmax
     &                         )
     
      call write_fit_time_series(file_fitted_out,itime,fitted,nobs,ndef,nloc)
      
      deallocate(obsdata)
      deallocate(fitted)
      deallocate(sig,ampci,phgi,ampc,phg,sdev,sdev0,rmsr,rmsr0,resmax,resmax_ls,iresmax,ssq,cov)
      
      
	end program versatile_analyze
	
	
      subroutine parse_command_line(buffer)
      character(LEN=*) buffer
      integer status
      call get_command_argument(1,buffer,status)
      buffer = trim(buffer)
      return
      end subroutine
	

