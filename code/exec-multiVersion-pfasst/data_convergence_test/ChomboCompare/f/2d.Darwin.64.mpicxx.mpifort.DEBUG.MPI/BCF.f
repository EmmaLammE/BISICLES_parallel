      subroutine EXTRAPGHOSTBC(
     &           state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,ibcBoxlo0,ibcBoxlo1
     &           ,ibcBoxhi0,ibcBoxhi1
     &           ,idir
     &           ,side
     &           ,dx
     &           ,startcomp
     &           ,endcomp
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL*8 state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer ibcBoxlo0,ibcBoxlo1
      integer ibcBoxhi0,ibcBoxhi1
      integer idir
      integer side
      REAL*8 dx
      integer startcomp
      integer endcomp
      integer nc
      integer ii0,ii1
      integer i0,i1
      REAL*8 nearval, farval
      ii0 = side*CHF_ID(0,idir)
      ii1 = side*CHF_ID(1,idir)
      do nc = startcomp, endcomp
      do i1 = ibcBoxlo1,ibcBoxhi1
      do i0 = ibcBoxlo0,ibcBoxhi0
         nearval = state(i0-ii0,i1-ii1,nc)
         farval  = state(i0-2*ii0,i1-2*ii1,nc)
         state(i0,i1,nc) = (2.0d0)*nearval - farval
      enddo
      enddo
      enddo
      return
      end
      subroutine HOEXTRAPGHOSTBC(
     &           state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,ibcBoxlo0,ibcBoxlo1
     &           ,ibcBoxhi0,ibcBoxhi1
     &           ,idir
     &           ,side
     &           ,dx
     &           ,startcomp
     &           ,endcomp
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL*8 state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer ibcBoxlo0,ibcBoxlo1
      integer ibcBoxhi0,ibcBoxhi1
      integer idir
      integer side
      REAL*8 dx
      integer startcomp
      integer endcomp
      integer nc
      integer ii0,ii1
      integer i0,i1
      REAL*8 nearval, midval, farval
      ii0 = side*CHF_ID(0,idir)
      ii1 = side*CHF_ID(1,idir)
      do nc = startcomp, endcomp
      do i1 = ibcBoxlo1,ibcBoxhi1
      do i0 = ibcBoxlo0,ibcBoxhi0
         nearval = state(i0-ii0,i1-ii1,nc)
         midval  = state(i0-2*ii0,i1-2*ii1,nc)
         farval  = state(i0-3*ii0,i1-3*ii1,nc)
         state(i0,i1,nc) = (3.0d0)*(nearval - midval) + farval
      enddo
      enddo
      enddo
      return
      end
      subroutine REFLECTGHOSTBC(
     &           state
     &           ,istatelo0,istatelo1
     &           ,istatehi0,istatehi1
     &           ,nstatecomp
     &           ,ibcBoxlo0,ibcBoxlo1
     &           ,ibcBoxhi0,ibcBoxhi1
     &           ,idir
     &           ,side
     &           ,dx
     &           ,startcomp
     &           ,endcomp
     &           ,scale
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL*8 state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           0:nstatecomp-1)
      integer ibcBoxlo0,ibcBoxlo1
      integer ibcBoxhi0,ibcBoxhi1
      integer idir
      integer side
      REAL*8 dx
      integer startcomp
      integer endcomp
      REAL*8 scale
      integer n
      integer ii0,ii1
      integer i0,i1
      if((side .ne. -1).and.(side.ne.1)) then
         call MAYDAYERROR()
      endif
      do n = startcomp, endcomp
      do i1 = ibcBoxlo1,ibcBoxhi1
      do i0 = ibcBoxlo0,ibcBoxhi0
         if (side.eq.-1)  then
            ii0 = CHF_ID(0,idir)*(2*(ibcBoxhi0-i0)+1)
            ii1 = CHF_ID(1,idir)*(2*(ibcBoxhi1-i1)+1)
         else if (side.eq.1) then
           ii0 = CHF_ID(0,idir)*(2*(ibcBoxlo0-i0)-1)
           ii1 = CHF_ID(1,idir)*(2*(ibcBoxlo1-i1)-1)
        else
           call MAYDAYERROR()
        endif
        state(i0,i1,n) = scale*state(i0+ii0,i1+ii1,n)
      enddo
      enddo
      enddo
      return
      end
