      subroutine ADJRHSSPEEDCTRL(
     & rhsx
     & ,irhsxlo0,irhsxlo1
     & ,irhsxhi0,irhsxhi1
     & ,rhsy
     & ,irhsylo0,irhsylo1
     & ,irhsyhi0,irhsyhi1
     & ,misfit
     & ,imisfitlo0,imisfitlo1
     & ,imisfithi0,imisfithi1
     & ,umx
     & ,iumxlo0,iumxlo1
     & ,iumxhi0,iumxhi1
     & ,umy
     & ,iumylo0,iumylo1
     & ,iumyhi0,iumyhi1
     & ,uox
     & ,iuoxlo0,iuoxlo1
     & ,iuoxhi0,iuoxhi1
     & ,uoy
     & ,iuoylo0,iuoylo1
     & ,iuoyhi0,iuoyhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irhsxlo0,irhsxlo1
      integer irhsxhi0,irhsxhi1
      REAL*8 rhsx(
     & irhsxlo0:irhsxhi0,
     & irhsxlo1:irhsxhi1)
      integer irhsylo0,irhsylo1
      integer irhsyhi0,irhsyhi1
      REAL*8 rhsy(
     & irhsylo0:irhsyhi0,
     & irhsylo1:irhsyhi1)
      integer imisfitlo0,imisfitlo1
      integer imisfithi0,imisfithi1
      REAL*8 misfit(
     & imisfitlo0:imisfithi0,
     & imisfitlo1:imisfithi1)
      integer iumxlo0,iumxlo1
      integer iumxhi0,iumxhi1
      REAL*8 umx(
     & iumxlo0:iumxhi0,
     & iumxlo1:iumxhi1)
      integer iumylo0,iumylo1
      integer iumyhi0,iumyhi1
      REAL*8 umy(
     & iumylo0:iumyhi0,
     & iumylo1:iumyhi1)
      integer iuoxlo0,iuoxlo1
      integer iuoxhi0,iuoxhi1
      REAL*8 uox(
     & iuoxlo0:iuoxhi0,
     & iuoxlo1:iuoxhi1)
      integer iuoylo0,iuoylo1
      integer iuoyhi0,iuoyhi1
      REAL*8 uoy(
     & iuoylo0:iuoyhi0,
     & iuoylo1:iuoyhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 uo,um,t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      uo = (uox(i0,i1)**(2.0d0) + uoy(i0,i1)**(2.0d0))**(0.500d0)
      um = (umx(i0,i1)**(2.0d0)+umy(i0,i1)**(2.0d0))**(0.500d0)
      misfit(i0,i1) = (uo-um)**(2.0d0)
      t = (uo/((1.0d-10) + um)-1.0)
      rhsx(i0,i1) = t * umx(i0,i1)
      rhsy(i0,i1) = t * umy(i0,i1)
      enddo
      enddo
      return
      end
      subroutine ADJRHSVELCTRL(
     & rhsx
     & ,irhsxlo0,irhsxlo1
     & ,irhsxhi0,irhsxhi1
     & ,rhsy
     & ,irhsylo0,irhsylo1
     & ,irhsyhi0,irhsyhi1
     & ,misfit
     & ,imisfitlo0,imisfitlo1
     & ,imisfithi0,imisfithi1
     & ,umx
     & ,iumxlo0,iumxlo1
     & ,iumxhi0,iumxhi1
     & ,umy
     & ,iumylo0,iumylo1
     & ,iumyhi0,iumyhi1
     & ,uox
     & ,iuoxlo0,iuoxlo1
     & ,iuoxhi0,iuoxhi1
     & ,uoy
     & ,iuoylo0,iuoylo1
     & ,iuoyhi0,iuoyhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irhsxlo0,irhsxlo1
      integer irhsxhi0,irhsxhi1
      REAL*8 rhsx(
     & irhsxlo0:irhsxhi0,
     & irhsxlo1:irhsxhi1)
      integer irhsylo0,irhsylo1
      integer irhsyhi0,irhsyhi1
      REAL*8 rhsy(
     & irhsylo0:irhsyhi0,
     & irhsylo1:irhsyhi1)
      integer imisfitlo0,imisfitlo1
      integer imisfithi0,imisfithi1
      REAL*8 misfit(
     & imisfitlo0:imisfithi0,
     & imisfitlo1:imisfithi1)
      integer iumxlo0,iumxlo1
      integer iumxhi0,iumxhi1
      REAL*8 umx(
     & iumxlo0:iumxhi0,
     & iumxlo1:iumxhi1)
      integer iumylo0,iumylo1
      integer iumyhi0,iumyhi1
      REAL*8 umy(
     & iumylo0:iumyhi0,
     & iumylo1:iumyhi1)
      integer iuoxlo0,iuoxlo1
      integer iuoxhi0,iuoxhi1
      REAL*8 uox(
     & iuoxlo0:iuoxhi0,
     & iuoxlo1:iuoxhi1)
      integer iuoylo0,iuoylo1
      integer iuoyhi0,iuoyhi1
      REAL*8 uoy(
     & iuoylo0:iuoyhi0,
     & iuoylo1:iuoyhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 uo,um,t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      rhsx(i0,i1) = uox(i0,i1) - umx(i0,i1)
      rhsy(i0,i1) = uoy(i0,i1) - umy(i0,i1)
      misfit(i0,i1) = (0.500d0) *
     & rhsx(i0,i1)**(2.0d0)
     & + rhsy(i0,i1)**(2.0d0)
      enddo
      enddo
      return
      end
      subroutine ADJRHSLOGSPDCTRL(
     & rhsx
     & ,irhsxlo0,irhsxlo1
     & ,irhsxhi0,irhsxhi1
     & ,rhsy
     & ,irhsylo0,irhsylo1
     & ,irhsyhi0,irhsyhi1
     & ,misfit
     & ,imisfitlo0,imisfitlo1
     & ,imisfithi0,imisfithi1
     & ,umx
     & ,iumxlo0,iumxlo1
     & ,iumxhi0,iumxhi1
     & ,umy
     & ,iumylo0,iumylo1
     & ,iumyhi0,iumyhi1
     & ,uox
     & ,iuoxlo0,iuoxlo1
     & ,iuoxhi0,iuoxhi1
     & ,uoy
     & ,iuoylo0,iuoylo1
     & ,iuoyhi0,iuoyhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irhsxlo0,irhsxlo1
      integer irhsxhi0,irhsxhi1
      REAL*8 rhsx(
     & irhsxlo0:irhsxhi0,
     & irhsxlo1:irhsxhi1)
      integer irhsylo0,irhsylo1
      integer irhsyhi0,irhsyhi1
      REAL*8 rhsy(
     & irhsylo0:irhsyhi0,
     & irhsylo1:irhsyhi1)
      integer imisfitlo0,imisfitlo1
      integer imisfithi0,imisfithi1
      REAL*8 misfit(
     & imisfitlo0:imisfithi0,
     & imisfitlo1:imisfithi1)
      integer iumxlo0,iumxlo1
      integer iumxhi0,iumxhi1
      REAL*8 umx(
     & iumxlo0:iumxhi0,
     & iumxlo1:iumxhi1)
      integer iumylo0,iumylo1
      integer iumyhi0,iumyhi1
      REAL*8 umy(
     & iumylo0:iumyhi0,
     & iumylo1:iumyhi1)
      integer iuoxlo0,iuoxlo1
      integer iuoxhi0,iuoxhi1
      REAL*8 uox(
     & iuoxlo0:iuoxhi0,
     & iuoxlo1:iuoxhi1)
      integer iuoylo0,iuoylo1
      integer iuoyhi0,iuoyhi1
      REAL*8 uoy(
     & iuoylo0:iuoyhi0,
     & iuoylo1:iuoyhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 uo,um,r,s,t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      uo = (uox(i0,i1)**(2.0d0)+uoy(i0,i1)**(2.0d0))**(0.500d0)
      um = (umx(i0,i1)**(2.0d0)+umy(i0,i1)**(2.0d0))**(0.500d0)
      uo = uo + 10.0
      um = um + 10.0
      r = um/uo
      t = log(r)
      s = t * uo / um**2
      rhsx(i0,i1) = s*uox(i0,i1)
      rhsy(i0,i1) = s*uoy(i0,i1)
      misfit(i0,i1) = (0.500d0) * t**2
      enddo
      enddo
      return
      end
      subroutine ADJRHSMASSCTRL(
     & rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,nrhscomp
     & ,misfit
     & ,imisfitlo0,imisfitlo1
     & ,imisfithi0,imisfithi1
     & ,thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,dx
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & 0:nrhscomp-1)
      integer imisfitlo0,imisfitlo1
      integer imisfithi0,imisfithi1
      REAL*8 misfit(
     & imisfitlo0:imisfithi0,
     & imisfitlo1:imisfithi1)
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      REAL*8 dx
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ip,jp
      integer dir,ndir
      REAL*8 oneontwodx
      oneontwodx = (0.500d0) / dx
      ndir = 2
      do dir =0,ndir-1
         ip = CHF_ID(0,dir)
         jp = CHF_ID(1,dir)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         rhs(i0,i1,dir) = thck(i0,i1) * oneontwodx
     & * (misfit(i0+ip,i1+jp)
     & - misfit(i0-ip,i1-jp))
      enddo
      enddo
      end do
      return
      end
      subroutine BOUNDEXPCTRL(
     & x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,a
     & ,ialo0,ialo1
     & ,iahi0,iahi1
     & ,lb
     & ,ub
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL*8 a(
     & ialo0:iahi0,
     & ialo1:iahi1)
      REAL*8 lb
      REAL*8 ub
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      t = min(ub,max(a(i0,i1),lb))
      x(i0,i1) = x(i0,i1) * exp(t)
      enddo
      enddo
      return
      end
      subroutine EXPCTRL(
     & x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,a
     & ,ialo0,ialo1
     & ,iahi0,iahi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL*8 a(
     & ialo0:iahi0,
     & ialo1:iahi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      x(i0,i1) = x(i0,i1)
     & * exp(a(i0,i1))
      enddo
      enddo
      return
      end
      subroutine BOUNDCTRL(
     & x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,lb
     & ,ub
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      REAL*8 lb
      REAL*8 ub
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      t = x(i0,i1)
      x(i0,i1) = min(ub,max(lb,t))
      enddo
      enddo
      return
      end
      subroutine INCRBOUNDCTRL(
     & z
     & ,izlo0,izlo1
     & ,izhi0,izhi1
     & ,y
     & ,iylo0,iylo1
     & ,iyhi0,iyhi1
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,lb
     & ,ub
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     & izlo0:izhi0,
     & izlo1:izhi1)
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     & iylo0:iyhi0,
     & iylo1:iyhi1)
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      REAL*8 lb
      REAL*8 ub
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      t = x(i0,i1)
      z(i0,i1) = y(i0,i1) + min(ub,max(lb,t))
      enddo
      enddo
      return
      end
      subroutine GRADBARRIERCTRL(
     & p
     & ,iplo0,iplo1
     & ,iphi0,iphi1
     & ,g
     & ,iglo0,iglo1
     & ,ighi0,ighi1
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,b
     & ,tol
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iplo0,iplo1
      integer iphi0,iphi1
      REAL*8 p(
     & iplo0:iphi0,
     & iplo1:iphi1)
      integer iglo0,iglo1
      integer ighi0,ighi1
      REAL*8 g(
     & iglo0:ighi0,
     & iglo1:ighi1)
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      REAL*8 b
      REAL*8 tol
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 xx,dm,dp,bsq
      bsq = b**2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      xx = x(i0,i1)
      if (xx.le.-b) then
         xx = -b+tol
      else if (xx.ge.b) then
         xx = b-tol
      end if
      dm = xx+b+tol
      dp = xx-b-tol
      p(i0,i1) = - log(- dp * dm / bsq)
      g(i0,i1) = - (1.0d0) / dm - (1.0d0) / dp
      enddo
      enddo
      return
      end
      subroutine MULTHATCTRL(
     & fab
     & ,ifablo0,ifablo1
     & ,ifabhi0,ifabhi1
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,lb
     & ,ub
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     & ifablo0:ifabhi0,
     & ifablo1:ifabhi1)
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      REAL*8 lb
      REAL*8 ub
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if ((x(i0,i1).le.lb).or.(x(i0,i1).ge.ub)) then
         fab(i0,i1) = (0.0d0)
      end if
      enddo
      enddo
      return
      end
      subroutine HARDPOINTINCTRL(
     & fab
     & ,ifablo0,ifablo1
     & ,ifabhi0,ifabhi1
     & ,x
     & ,ixlo0,ixlo1
     & ,ixhi0,ixhi1
     & ,lb
     & ,ub
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     & ifablo0:ifabhi0,
     & ifablo1:ifabhi1)
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     & ixlo0:ixhi0,
     & ixlo1:ixhi1)
      REAL*8 lb
      REAL*8 ub
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (x(i0,i1).le.lb) then
         fab(i0,i1) = min((0.0d0),fab(i0,i1))
      else if (x(i0,i1).ge.ub) then
         fab(i0,i1) = max((0.0d0),fab(i0,i1))
      end if
      enddo
      enddo
      return
      end
      subroutine CADOTBCTRL(
     & r
     & ,irlo0,irlo1
     & ,irhi0,irhi1
     & ,c
     & ,iclo0,iclo1
     & ,ichi0,ichi1
     & ,a
     & ,ialo0,ialo1
     & ,iahi0,iahi1
     & ,nacomp
     & ,b
     & ,iblo0,iblo1
     & ,ibhi0,ibhi1
     & ,nbcomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer irlo0,irlo1
      integer irhi0,irhi1
      REAL*8 r(
     & irlo0:irhi0,
     & irlo1:irhi1)
      integer iclo0,iclo1
      integer ichi0,ichi1
      REAL*8 c(
     & iclo0:ichi0,
     & iclo1:ichi1)
      integer nacomp
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL*8 a(
     & ialo0:iahi0,
     & ialo1:iahi1,
     & 0:nacomp-1)
      integer nbcomp
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      REAL*8 b(
     & iblo0:ibhi0,
     & iblo1:ibhi1,
     & 0:nbcomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, icomp
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      r(i0,i1) = c(i0,i1) *
     & a(i0,i1,0) *
     & b(i0,i1,0)
      enddo
      enddo
      if (nacomp.gt.1) then
         do icomp = 1, nacomp - 1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            r(i0,i1) = r(i0,i1) +
     & c(i0,i1) *
     & a(i0,i1,icomp) *
     & b(i0,i1,icomp)
      enddo
      enddo
         end do
      end if
      return
      end
      subroutine LIMITINCRCTRL(
     & h
     & ,ihlo0,ihlo1
     & ,ihhi0,ihhi1
     & ,dh
     & ,idhlo0,idhlo1
     & ,idhhi0,idhhi1
     & ,ho
     & ,iholo0,iholo1
     & ,ihohi0,ihohi1
     & ,limit
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     & ihlo0:ihhi0,
     & ihlo1:ihhi1)
      integer idhlo0,idhlo1
      integer idhhi0,idhhi1
      REAL*8 dh(
     & idhlo0:idhhi0,
     & idhlo1:idhhi1)
      integer iholo0,iholo1
      integer ihohi0,ihohi1
      REAL*8 ho(
     & iholo0:ihohi0,
     & iholo1:ihohi1)
      REAL*8 limit
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      t = h(i0,i1) + dh(i0,i1)
      t = max(t,(0.0d0))
      t = max(ho(i0,i1)-limit,t)
      t = min(ho(i0,i1)+limit,t)
      h(i0,i1) = t
      enddo
      enddo
      return
      end
      subroutine UPDATEHCTRLB(
     & h
     & ,ihlo0,ihlo1
     & ,ihhi0,ihhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,dh
     & ,idhlo0,idhlo1
     & ,idhhi0,idhhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,ho
     & ,iholo0,iholo1
     & ,ihohi0,ihohi1
     & ,r
     & ,dhmax
     & ,hmin
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     & ihlo0:ihhi0,
     & ihlo1:ihhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      integer idhlo0,idhlo1
      integer idhhi0,idhhi1
      REAL*8 dh(
     & idhlo0:idhhi0,
     & idhlo1:idhhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer iholo0,iholo1
      integer ihohi0,ihohi1
      REAL*8 ho(
     & iholo0:ihohi0,
     & iholo1:ihohi1)
      REAL*8 r
      REAL*8 dhmax
      REAL*8 hmin
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 d, t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      dh(i0,i1) = max(dh(i0,i1),-h(i0,i1))
      if ( h(i0,i1) .lt. (1.0d-2)) then
         dh(i0,i1) = (0.0d0)
      end if
      t = h(i0,i1)
      h(i0,i1) = h(i0,i1) + dh(i0,i1)
      if ( h(i0,i1) .lt. hmin) then
         h(i0,i1) = (0.0d0)
      end if
      if ( h(i0,i1) .gt. ho(i0,i1) + dhmax) then
         h(i0,i1) = ho(i0,i1) + dhmax
      else if (h(i0,i1) .lt. ho(i0,i1) - dhmax) then
         h(i0,i1) = ho(i0,i1) - dhmax
      end if
      dh(i0,i1) = h(i0,i1) - t
      if (mask (i0,i1) .eq. (2)) then
         d = min(-dh(i0,i1), (0.0d0))
      else if (mask (i0,i1) .eq. (1)) then
         d = -dh(i0,i1)
      else
         d = (0.0d0)
      end if
      topg(i0,i1) = topg(i0,i1) + d
      enddo
      enddo
      return
      end
      subroutine UPDATEHCTRLC(
     & h
     & ,ihlo0,ihlo1
     & ,ihhi0,ihhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,dh
     & ,idhlo0,idhlo1
     & ,idhhi0,idhhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,ho
     & ,iholo0,iholo1
     & ,ihohi0,ihohi1
     & ,r
     & ,dhmax
     & ,hmin
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     & ihlo0:ihhi0,
     & ihlo1:ihhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      integer idhlo0,idhlo1
      integer idhhi0,idhhi1
      REAL*8 dh(
     & idhlo0:idhhi0,
     & idhlo1:idhhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer iholo0,iholo1
      integer ihohi0,ihohi1
      REAL*8 ho(
     & iholo0:ihohi0,
     & iholo1:ihohi1)
      REAL*8 r
      REAL*8 dhmax
      REAL*8 hmin
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 d, t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      dh(i0,i1) = max(dh(i0,i1),-h(i0,i1))
      if ( h(i0,i1) .lt. (1.0d-2)) then
         dh(i0,i1) = (0.0d0)
      end if
      t = h(i0,i1)
      h(i0,i1) = h(i0,i1) + dh(i0,i1)
      if ( h(i0,i1) .lt. hmin) then
         h(i0,i1) = (0.0d0)
      end if
      if ( h(i0,i1) .gt. ho(i0,i1) + dhmax) then
         h(i0,i1) = ho(i0,i1) + dhmax
      else if (h(i0,i1) .lt. ho(i0,i1) - dhmax) then
         h(i0,i1) = ho(i0,i1) - dhmax
      end if
      dh(i0,i1) = h(i0,i1) - t
      if (mask (i0,i1) .eq. (2)) then
         d = min(-dh(i0,i1), (0.0d0))
      else if (mask (i0,i1) .eq. (1)) then
         d = min(-dh(i0,i1), (0.0d0))
      else
         d = (0.0d0)
      end if
      topg(i0,i1) = topg(i0,i1) + d
      enddo
      enddo
      return
      end
      subroutine CONVOLVECTRLRB(
     & u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,nucomp
     & ,rb
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1,
     & 0:nucomp-1)
      integer rb
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      integer n,ncomp,imin,imax,indtot
      REAL*8 w,t,twow
      w = (2.0d0) * 2;
      twow = (2.0d0) * w;
      ncomp = nucomp
      do n = 0, ncomp - 1
            do j=iboxlo1, iboxhi1
               imin = iboxlo0
               indtot = imin + j
               imin = imin + abs(mod(indtot + rb, 2))
               imax = iboxhi0
               do i = imin, imax, 2
                  t = (
     & u(i+1,j,n)
     & + u(i-1,j,n)
     & + u(i,j+1,n)
     & + u(i,j-1,n)
     & )
                  u(i,j,n) = (t + w *u(i,j,n))/(twow)
               end do
            end do
      end do
      return
      end
      subroutine CONVOLVECTRL(
     & u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,nucomp
     & ,r
     & ,irlo0,irlo1
     & ,irhi0,irhi1
     & ,nrcomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1,
     & 0:nucomp-1)
      integer nrcomp
      integer irlo0,irlo1
      integer irhi0,irhi1
      REAL*8 r(
     & irlo0:irhi0,
     & irlo1:irhi1,
     & 0:nrcomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      integer n,ncomp,imin,imax,indtot
      REAL*8 w,t,twow
      w = (2.0d0) * 2;
      twow = (2.0d0) * w;
      ncomp = nucomp
      do n = 0, ncomp - 1
            do j=iboxlo1, iboxhi1
               do i=iboxlo0, iboxhi0
                  t = (
     & r(i+1,j,n)
     & + r(i-1,j,n)
     & + r(i,j+1,n)
     & + r(i,j-1,n)
     & )
                  u(i,j,n) = (t + w*r(i,j,n))/(twow)
               end do
            end do
      end do
      return
      end
