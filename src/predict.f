      
      subroutine vt_predict(itime,
     1                     ntime,
     1                     values,
     1                     ndef,
     1                     nloc,
     1                     latd, latm, 
     1                     istart_ha,iend_ha,
     1                     nconstituent,
     1                     itrend,
     1                     all_names,
     1                     all_freq,
     1                     all_amp,
     1                     all_phase)
     
      use constituent
      implicit none
      integer, parameter :: max_cnstnt = 80
      ! args
      integer,intent(in )::  istart_ha(7)
      integer,intent(in) :: iend_ha(7)
      integer,intent(in) :: itime(7,ntime) !< times at which tides are to be produced
      integer,intent(in) :: ntime          !< number of times at which output is requested
      real*8,intent(out) :: values(ntime,ndef,nloc) !< output values of prediction    
      integer,intent(in) :: ndef           !< number of components, 1=scalar, 2=vector (e.g. velocity)
      integer,intent(in) :: nloc           !< number of locations or stations
      integer,intent(in) :: latd           !< representative latitude (degree part)
      integer,intent(in) :: latm           !< representative longitude (minute part)
      integer,intent(in)  :: nconstituent   !< number of constituents (main and inferred are not distinguished)
      integer,intent(in)  :: itrend         !< whether or not a trend is included, 0=no and 1=yes !todo: trend is not accommodated

      character*5, dimension(max_cnstnt) :: all_names !< names of all constituents
      real*8,intent(in) :: all_freq(max_cnstnt)                  !< frequencies of all constituents (not used)
      real*8,intent(in) :: all_amp(max_cnstnt,ndef,nloc)         !< amplitudes of all constituents
      real*8,intent(in) :: all_phase(max_cnstnt,ndef,nloc)       !< phase of all constituents

      ! locals
      real*8 :: all_phase_rad(max_cnstnt,ndef,nloc)
      real*8 time(ntime)
      integer, parameter :: lcf = 58
      integer :: i,j,k
      real*8 xlat
      integer ic1,iy1,im1,id1,ic2,iy2,im2,id2
      integer kd(ntime)
      integer kd_anal1   ! day at beginning of original analysis
      integer kd_anal2   ! day at end of original analysis
      real*8 x(ntime)
      integer kd1,kd2,kh1,kh2,khm
      real*8 cenhr,cumhr,xmid
      real*8 anal_hr1, anal_hr2
      integer isec  ! todo???
      integer kdtp(ntime) !todo: what is this?
      real*8 rel_time(ntime) !todo: not same time as it was in analysis
      integer :: kh
      real*8 :: hr
      real*8 :: sum_components(ndef,nloc)
      real*8 :: vx,ux,fx,arg(ndef,nloc)
      
      ! decimal location
      xlat=latd+latm/60.
	
      
      ! convert to gregorian day
      iy1 = itime(1,1)
      im1 = itime(2,1)
      id1 = itime(3,1)
      ic1 = itime(7,1)
      iy2 = itime(1,ntime)
      im2 = itime(2,ntime)
      id2 = itime(3,ntime)
      ic2 = itime(7,ntime)
      do i = 1,ntime
        call GDAY(itime(3,i),itime(2,i),itime(1,i),itime(7,i),kd(i))
      end do
      
      
      
      ! Calculate day and hour for start/end/midpoint of analysis time period
      ! this is the time centering used for the trend. 
      ! Probably should make this a single subroutine so that it can't get out of sync
      call GDAY(istart_ha(3),istart_ha(2),istart_ha(1),istart_ha(7),kd_anal1)
      call GDAY(iend_ha(3),iend_ha(2),iend_ha(1),iend_ha(7),kd_anal2)

      kd1=kd_anal1
      kd2=kd_anal2
      kh1=24*kd1
      kh2=24*(kd2+1)
      khm=(kh1+kh2)/2    !hour at midpoint of original analysis

      cenhr=dfloat((kh2-kh1)/2) ! division by two is exact because these are multiples of 24
                                ! cenhr may not be the dead center of the analysis in edge cases
      anal_hr1=-cenhr+24.d0*0.d0+istart_ha(5)/60.d0+istart_ha(6)/3600.d0
      anal_hr2=-cenhr+24.d0*(kd_anal2-kd_anal1)+iend_ha(5)/60.d0+iend_ha(6)/3600.d0
      xmid = 5.d-1*(anal_hr1+anal_hr2)
      

      all_phase_rad = all_phase/360.d0
      
	k=1
	isec=0
	do 20 i=1,ntime
        cumhr=-cenhr+24.D0*(kd(i)-kd1)
        time(i)=cumhr+dfloat(itime(4,i))+(dfloat(itime(5,i))
     1             +dfloat(isec)/60.d0)/60.d0
        x(k)=time(i)-anal_hr1
        kdtp(k)=kd(i)    ! todo: see if two variables are needed
        k=k+1
20    continue      

      ! todo: time(i) seems to already be redundant with rel_time(i). 
      ! The code seems to convert it to x (time since first timestamp) only
      ! to convert it back to rel_time which is with respect to the center
      
      do i=1,ntime
        rel_time(i) = (x(i)-xmid)/(24.*365.)  ! for use with itrend, was in q()
	  kh=24*kdtp(i)+itime(4,i)
	  hr=dfloat(kh)+itime(5,i)/60.d0+itime(6,i)/3600.d0
        call SETVUF(hr,xlat,nconstituent,all_names)
        sum_components=all_amp(1,:,:)
        if(itrend /= 0) then
            ! note that the coefficient for the trend is stored as the phase of Z0
            ! todo: test this gives the same as original data
             sum_components=sum_components+all_phase(1,:,:)*(x(i)-xmid)/(365.d0*24.d0)
        end if
        do j=2,nconstituent
          call VUF(all_names(j),vx,ux,fx)
          arg=(vx+ux-all_phase_rad(j,:,:))*twopi
          sum_components=sum_components+fx*all_amp(j,:,:)*cos(arg)
        end do  ! j
        values(i,:,:) = sum_components
      end do  ! i
      
      return
      end subroutine