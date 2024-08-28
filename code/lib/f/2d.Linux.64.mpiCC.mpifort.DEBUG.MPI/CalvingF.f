      subroutine EFFECTIVETHICKNESS(
     & he
     & ,ihelo0,ihelo1
     & ,ihehi0,ihehi1
     & ,h
     & ,ihlo0,ihlo1
     & ,ihhi0,ihhi1
     & ,m
     & ,imlo0,imlo1
     & ,imhi0,imhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ihelo0,ihelo1
      integer ihehi0,ihehi1
      REAL*8 he(
     & ihelo0:ihehi0,
     & ihelo1:ihehi1)
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     & ihlo0:ihhi0,
     & ihlo1:ihhi1)
      integer imlo0,imlo1
      integer imhi0,imhi1
      REAL*8 m(
     & imlo0:imhi0,
     & imlo1:imhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 mm
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      mm = min(m(i0,i1),(1.0d0))
      if ( mm .gt. (0.0d0) ) then
         he(i0,i1) = h(i0,i1)/ mm
      else
         he(i0,i1) = (0.0d0)
      end if
      enddo
      enddo
      return
      end
