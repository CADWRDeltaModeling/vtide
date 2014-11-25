
!****************************************************************************
!
!  PROGRAM: versatile_predict
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
      integer :: nobs
      integer, parameter ::max_obs = 106000
      integer,allocatable :: itime(:,:)
      integer :: ipredtime(7,max_obs)
      integer :: istart_ha(7)
      integer :: iend_ha(7)
      real*8, allocatable :: obsdata(:,:,:)

!c    For prediction
      character*80 :: fname = "all_constituent.out"
      integer :: station_no
      integer kd1,kd2
      character*20 :: station_name_pred
      integer :: latd,latm
      integer :: lond,lonm
      integer :: nconstituent
      integer :: ntvl_min
      character*5, dimension(max_main_cnstnt) :: all_names !todo: should be all constituents
      real*8, dimension(max_main_cnstnt,2,3) :: all_freq
      real*8, dimension(max_main_cnstnt,2,3) :: all_amp
      real*8, dimension(max_main_cnstnt,2,3) :: all_phase
      real*8,allocatable :: values(:,:,:)
      integer i
      integer :: idd,imm,iyy,icc,ikd
      integer nperday,istep,imin,ihh,isec,ntime
      call parse_command_line(main_file)
     
      open(UNIT=KIN,FILE=main_file,STATUS='OLD',ACTION='READ')
      read(KIN,'(a)') file_const_db       !6
      read(KIN,'(a)') file_analysis
      read(KIN,'(4I2)') im1,id1,ic1,iy1
      read(KIN,'(4I2)') im2,id2,ic2,iy2
      read(KIN,*),ndef,nloc
      read(KIN,*) ntvl_min
      close(KIN)
      
      call GDAY(ID1,IM1,IY1,IC1,kd1)
      call GDAY(ID2,IM2,IY2,IC2,kd2)
      nperday=(24*60)/ntvl_min
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
      file_obs='kiw06.dat'
      nobs=max_obs
      if (.not. allocated(obsdata)) allocate(obsdata(nobs,ndef,nloc))
      if (.not. allocated(itime)) allocate(itime(7,nobs))   
      call vt_read_constituent_db(file_const_db)      
      call vt_read_data(file_obs,ndef,nloc,
     &                  id1,im1,iy1,ic1,
     &                  id2,im2,iy2,ic2,
     &                  nobs,itime,obsdata)
      
      
      call vt_read_constituent_file(fname, 
     1                         station_no, 
     1                         station_name_pred, 
     1                         latd, latm, 
     1                         lond,lonm,
     1                         istart_ha,iend_ha,
     1                         nconstituent,ndef,nloc,itrend,
     1                         all_names,
     1                         all_freq,
     1                         all_amp,
     1                         all_phase)
     
      
     
      if (.not. allocated(values)) allocate(values(nobs,ndef,nloc))
      
      call vt_predict(ipredtime,
     1                ntime,
     1                values,
     1                ndef,
     1                nloc,
     1                latd, latm, 
     1                istart_ha,iend_ha,
     1                nconstituent,
     1                itrend,
     1                all_names,
     1                all_freq,
     1                all_amp,
     1                all_phase)
      

      
      deallocate(obsdata)
	end program versatile_predict
	
	
      subroutine parse_command_line(buffer)
      character(LEN=*) buffer
      integer status
      call get_command_argument(1,buffer,status)
      buffer = trim(buffer)
      return
      end subroutine
	

