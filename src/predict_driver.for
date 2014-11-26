
!****************************************************************************
!
!  PROGRAM: web
!
!  PURPOSE:  Main shell program.
!
!****************************************************************************

      program versatile_predict
      use file_io
      use constituent
	  implicit none
	  integer , parameter :: KIN = 4
      character(LEN=80) :: main_file
      integer, parameter :: max_main_cnstnt = 80
      integer, parameter :: max_infer = 80
      integer :: ndef
      integer :: nloc
      integer :: itrend
      integer :: id1,im1,iy1,id2,im2,iy2,ic1,ic2

      integer, allocatable :: ipredtime(:,:)      

!c    For prediction
      character*80 :: file_predict 
      integer :: station_no
      integer kd1,kd2
      character*20 :: station_name_pred
      integer :: latd,latm
      integer :: lond,lonm
      integer :: nconstituent
      integer :: ntvl_min
      character*5, dimension(max_main_cnstnt) :: all_names !todo: should be all constituents
      real*8, allocatable,dimension(:) :: all_freq
      real*8, allocatable,dimension(:,:,:) :: all_amp
      real*8, allocatable,dimension(:,:,:) :: all_phase
      real*8,allocatable :: values(:,:,:)
      integer i
      integer :: idd,imm,iyy,icc,ikd
      integer nperday,istep,imin,ihh,isec,ntime
      character :: sep1,sep2
      integer ::  istart_ha(7)  ! timestamp of original analysis start (for trend)
      integer ::  iend_ha(7)    ! timestamp of original analysis end
      
      
      
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
!c	Original, fitted, and residual time series are output to file 25 while
!c	the same are also output to file 26 in a format that could be input to
!c	Excel for plotting.
!c
      call parse_command_line(main_file)
     
      open(UNIT=KIN,FILE=main_file,STATUS='OLD',ACTION='READ')
      read(KIN,*) file_constituent_out
      read(KIN,*) file_const_db
      read(KIN,*) file_predict
      read(KIN,'(I2,I2,1a,I2,1a,I2)') ic1,iy1,sep1,im1,sep2,id1
      read(KIN,'(I2,I2,1a,I2,1a,I2)') ic2,iy2,sep1,im2,sep2,id2
      read(KIN,*) ntvl_min
      close(KIN)
      
      call GDAY(ID1,IM1,IY1,IC1,kd1)
      call GDAY(ID2,IM2,IY2,IC2,kd2)
      nperday=(24*60)/ntvl_min          ! todo: mention noneven intervals
      ntime = ((kd2-kd1) + 1) * nperday ! todo: is this calculation correct?
      if (.not. allocated(ipredtime)) allocate(ipredtime(7,ntime))
      
      i=0
      do ikd=kd1,kd2
          call dmy(IDD,IMM,IYY,ICC,ikd) 
          do istep=1,nperday
              i=i+1
              imin=(istep-1)*ntvl_min
              ihh=imin/60
              imin=mod(imin,60)
              isec = 0
              ipredtime(1,i)=iyy
              ipredtime(2,i)=imm
              ipredtime(3,i)=idd
              ipredtime(4,i)=ihh
              ipredtime(5,i)=imin
              ipredtime(6,i)=isec
	        ipredtime(7,i)=icc
          end do
      end do
      ntime=i

      call vt_read_constituent_db(file_const_db) 

      ! At the moment we are really only calling this to get ndef and nloc for allocation
      ! todo: could make the number of constituents dynamic as well.
      call vt_read_constituent_file_header(file_constituent_out, 
     1                         station_no, 
     1                         station_name_pred, 
     1                         latd, latm, 
     1                         lond,lonm,
     1                         istart_ha,iend_ha,
     1                         nconstituent,ndef,nloc,itrend) 

      if (.not. allocated(all_freq)) allocate(all_freq(max_main_cnstnt))
      if (.not. allocated(all_amp)) allocate(all_amp(max_main_cnstnt,ndef,nloc))
      if (.not. allocated(all_phase)) allocate(all_phase(max_main_cnstnt,ndef,nloc))
      
      call vt_read_constituent_file(file_constituent_out, 
     1                         station_no, 
     1                         station_name_pred, 
     1                         latd, latm, 
     1                         lond,lonm,
     1                         istart_ha,
     1                         iend_ha,  
     1                         nconstituent,ndef,nloc,itrend,
     1                         all_names,
     1                         all_freq,
     1                         all_amp,
     1                         all_phase)    
     
      
     
      if (.not. allocated(values)) allocate(values(ntime,ndef,nloc))
      
      call vt_predict(ipredtime,
     1                ntime,
     1                values,
     1                ndef,
     1                nloc,
     1                latd, latm, 
     1                istart_ha,
     1                iend_ha,
     1                nconstituent,
     1                itrend,
     1                all_names,
     1                all_freq,
     1                all_amp,
     1                all_phase)
      
      call write_fit_time_series(file_predict,ipredtime,values,ntime,ndef,nloc)

      deallocate(all_freq)
      deallocate(all_phase)
      deallocate(all_amp)
 !c     deallocate(obsdata)
	end program versatile_predict
	
	
      subroutine parse_command_line(buffer)
      character(LEN=*) buffer
      integer status
      call get_command_argument(1,buffer,status)
      buffer = trim(buffer)
      return
      end subroutine
	

