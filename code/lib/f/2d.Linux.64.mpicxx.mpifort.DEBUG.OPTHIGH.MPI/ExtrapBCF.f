      subroutine SIMPLEEXTRAPBC(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,ibcboxlo0,ibcboxlo1
     & ,ibcboxhi0,ibcboxhi1
     & ,dir
     & ,hiLo
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer ibcboxlo0,ibcboxlo1
      integer ibcboxhi0,ibcboxhi1
      integer dir
      integer hiLo
      integer i0,i1
      integer ii0,ii1
      integer n
      integer offset
      offset = 1
      if (hiLo.eq.0) offset = -1
      ii0 = offset*CHF_ID(0,dir)
                ii1 = offset*CHF_ID(1,dir)
      do n = 0, nphicomp-1
      do i1 = ibcboxlo1,ibcboxhi1
      do i0 = ibcboxlo0,ibcboxhi0
            phi(i0,i1,n) = (2.0d0)*phi(i0-ii0,i1-ii1,n)
     & - phi(i0-2*ii0,i1-2*ii1,n)
      enddo
      enddo
      enddo
      return
      end
      subroutine OLDSIMPLEREFLECTBC(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,ibcboxlo0,ibcboxlo1
     & ,ibcboxhi0,ibcboxhi1
     & ,dir
     & ,hiLo
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer ibcboxlo0,ibcboxlo1
      integer ibcboxhi0,ibcboxhi1
      integer dir
      integer hiLo
      integer i0,i1
      integer ii0,ii1
      integer n
      integer offset
      offset = 1
      if (hiLo.eq.0) offset = -1
      ii0 = offset*CHF_ID(0,dir)
                ii1 = offset*CHF_ID(1,dir)
      do n = 0, nphicomp-1
      do i1 = ibcboxlo1,ibcboxhi1
      do i0 = ibcboxlo0,ibcboxhi0
            phi(i0,i1,n) = phi(i0-2*ii0,i1-2*ii1,n)
      enddo
      enddo
      enddo
      return
      end
      subroutine SIMPLEREFLECTBC(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,ibcboxlo0,ibcboxlo1
     & ,ibcboxhi0,ibcboxhi1
     & ,dir
     & ,hiLo
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer ibcboxlo0,ibcboxlo1
      integer ibcboxhi0,ibcboxhi1
      integer dir
      integer hiLo
      integer i0,i1
      integer ii0,ii1
      integer n
      integer offset
      offset = 1
      if (hiLo.eq.0) offset = -1
      ii0 = offset*CHF_ID(0,dir)
                ii1 = offset*CHF_ID(1,dir)
      do n = 0, nphicomp-1
      do i1 = ibcboxlo1,ibcboxhi1
      do i0 = ibcboxlo0,ibcboxhi0
            phi(i0,i1,n) = phi(i0-1*ii0,i1-1*ii1,n)
      enddo
      enddo
      enddo
      return
      end
      subroutine EXTRAPCORNER2D(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer n
      do n = 0, nphicomp -1
         phi(iboxhi0+1,iboxhi1+1, n) = - phi(iboxhi0,iboxhi1, n)
     & + phi(iboxhi0,iboxhi1+1,n) + phi(iboxhi0+1,iboxhi1,n)
         phi(iboxhi0+1,iboxlo1-1, n) = - phi(iboxhi0,iboxlo1, n)
     & + phi(iboxhi0,iboxlo1-1,n) + phi(iboxhi0+1,iboxlo1,n)
         phi(iboxlo0-1,iboxlo1-1, n) = - phi(iboxlo0,iboxlo1, n)
     & + phi(iboxlo0-1,iboxlo1,n) + phi(iboxlo0,iboxlo1-1,n)
         phi(iboxlo0-1,iboxhi1+1, n) = - phi(iboxlo0,iboxhi1, n)
     & + phi(iboxlo0,iboxhi1+1,n) + phi(iboxlo0-1,iboxhi1,n)
      end do
      return
      end
