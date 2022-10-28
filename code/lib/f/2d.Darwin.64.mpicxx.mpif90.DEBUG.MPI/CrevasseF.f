      subroutine VDVSTRESSB(
     &           k
     &           ,iklo0,iklo1
     &           ,ikhi0,ikhi1
     &           ,h
     &           ,ihlo0,ihlo1
     &           ,ihhi0,ihhi1
     &           ,hp
     &           ,ihplo0,ihplo1
     &           ,ihphi0,ihphi1
     &           ,d
     &           ,idlo0,idlo1
     &           ,idhi0,idhi1
     &           ,rxx
     &           ,irxxlo0,irxxlo1
     &           ,irxxhi0,irxxhi1
     &           ,rhoi
     &           ,rhow
     &           ,gravity
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklo0,iklo1
      integer ikhi0,ikhi1
      REAL*8 k(
     &           iklo0:ikhi0,
     &           iklo1:ikhi1)
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     &           ihlo0:ihhi0,
     &           ihlo1:ihhi1)
      integer ihplo0,ihplo1
      integer ihphi0,ihphi1
      REAL*8 hp(
     &           ihplo0:ihphi0,
     &           ihplo1:ihphi1)
      integer idlo0,idlo1
      integer idhi0,idhi1
      REAL*8 d(
     &           idlo0:idhi0,
     &           idlo1:idhi1)
      integer irxxlo0,irxxlo1
      integer irxxhi0,irxxhi1
      REAL*8 rxx(
     &           irxxlo0:irxxhi0,
     &           irxxlo1:irxxhi1)
      REAL*8 rhoi
      REAL*8 rhow
      REAL*8 gravity
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 sqrtpid, pi, a, b, c, ires,ierr, ki, kw
      pi = 4.0 * atan(1.0) 
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (h(i0,i1).gt.(0.0d0)) then
         sqrtpid = (d(i0,i1) * pi)**(0.500d0)
         a = rxx(i0,i1)-rhoi*gravity*h(i0,i1)
         b = rhoi*gravity 
         call vdvint(ires,ierr,(0.0d0), d(i0,i1),
     &        d(i0,i1),h(i0,i1),a,b)
         k(i0,i1) = (2.0d0) / sqrtpid * ires
         if (hp(i0,i1) .gt. (0.0d0)) then
            c = min(hp(i0,i1),d(i0,i1))
            a = rhow*gravity*hp(i0,i1)
            b = -rhow*gravity 
            call vdvint(ires,ierr,(0.0d0),c,
     &           d(i0,i1),h(i0,i1),a,b)
            kw = (2.0d0) / sqrtpid * ires
            k(i0,i1) = k(i0,i1) + kw
         end if
      else
         k(i0,i1) = (0.0d0)
      end if
      enddo
      enddo
      return
      end
      subroutine VDVSTRESSS(
     &           k
     &           ,iklo0,iklo1
     &           ,ikhi0,ikhi1
     &           ,h
     &           ,ihlo0,ihlo1
     &           ,ihhi0,ihhi1
     &           ,dw
     &           ,idwlo0,idwlo1
     &           ,idwhi0,idwhi1
     &           ,d
     &           ,idlo0,idlo1
     &           ,idhi0,idhi1
     &           ,rxx
     &           ,irxxlo0,irxxlo1
     &           ,irxxhi0,irxxhi1
     &           ,rhoi
     &           ,rhow
     &           ,gravity
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklo0,iklo1
      integer ikhi0,ikhi1
      REAL*8 k(
     &           iklo0:ikhi0,
     &           iklo1:ikhi1)
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     &           ihlo0:ihhi0,
     &           ihlo1:ihhi1)
      integer idwlo0,idwlo1
      integer idwhi0,idwhi1
      REAL*8 dw(
     &           idwlo0:idwhi0,
     &           idwlo1:idwhi1)
      integer idlo0,idlo1
      integer idhi0,idhi1
      REAL*8 d(
     &           idlo0:idhi0,
     &           idlo1:idhi1)
      integer irxxlo0,irxxlo1
      integer irxxhi0,irxxhi1
      REAL*8 rxx(
     &           irxxlo0:irxxhi0,
     &           irxxlo1:irxxhi1)
      REAL*8 rhoi
      REAL*8 rhow
      REAL*8 gravity
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 sqrtpid, pi, a, b, c, ires,ierr, ki,kw
      pi = 4.0 * atan(1.0) 
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (h(i0,i1).gt.(0.0d0)) then
         sqrtpid = (d(i0,i1) * pi)**(0.500d0)
         a = rxx(i0,i1)
         b = -rhoi*gravity 
         call vdvint(ires,ierr,(0.0d0), d(i0,i1),
     &        d(i0,i1),h(i0,i1),a,b)
         k(i0,i1) = (2.0d0) / sqrtpid * ires
         if (dw(i0,i1) .gt. (0.0d0)) then
            c = d(i0,i1)-dw(i0,i1)
            a = -rhoi*gravity*c 
            b = rhow*gravity 
            call vdvint(ires,ierr,c,d(i0,i1),
     &           d(i0,i1),h(i0,i1),a,b)
            kw = (2.0d0) / sqrtpid * ires
            k(i0,i1) = k(i0,i1) + kw
         end if
      else
         k(i0,i1) = (0.0d0)
      end if
      enddo
      enddo
      return
      end
