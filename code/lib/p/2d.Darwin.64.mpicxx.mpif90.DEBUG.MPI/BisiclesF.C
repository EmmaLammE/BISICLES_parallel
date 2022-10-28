#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine ENFORCEWELLPOSEDCELL(
     &           C
     &           ,iClo0,iClo1
     &           ,iChi0,iChi1
     &           ,mux
     &           ,imuxlo0,imuxlo1
     &           ,imuxhi0,imuxhi1
     &           ,muy
     &           ,imuylo0,imuylo1
     &           ,imuyhi0,imuyhi1
     &           ,mu0
     &           ,C0
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iClo0,iClo1
      integer iChi0,iChi1
      REAL_T C(
     &           iClo0:iChi0,
     &           iClo1:iChi1)
      integer imuxlo0,imuxlo1
      integer imuxhi0,imuxhi1
      REAL_T mux(
     &           imuxlo0:imuxhi0,
     &           imuxlo1:imuxhi1)
      integer imuylo0,imuylo1
      integer imuyhi0,imuyhi1
      REAL_T muy(
     &           imuylo0:imuyhi0,
     &           imuylo1:imuyhi1)
      REAL_T mu0
      REAL_T C0
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
       integer i0,i1
       integer ix0,ix1
       integer iy0,iy1
       Real_T musum
       ix0 = CHF_ID(0,0)
                 ix1 = CHF_ID(0,1)
      iy0 = CHF_ID(1,0)
                iy1 = CHF_ID(1,1)          
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      if (C(i0,i1).lt.C0) then
         musum = mux(i0,i1) +  mux(i0+ix0,i1+ix1) 
     &        +  muy(i0,i1) +  muy(i0+iy0,i1+iy1)
         if (musum.lt.mu0) then
            C(i0,i1) = C0
         end if
      end if
      
      enddo
      enddo
      return
      end
      subroutine MASKEDREPLACE(
     &           a
     &           ,ialo0,ialo1
     &           ,iahi0,iahi1
     &           ,b
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,m
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1)
      integer m
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         if (mask(i0,i1).eq.m) then
            a(i0,i1) =  b(i0,i1)
         end if
      
      enddo
      enddo
      return
      end
      subroutine PROXIMITYFILL(
     &           a
     &           ,ialo0,ialo1
     &           ,iahi0,iahi1
     &           ,b
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,p
     &           ,iplo0,iplo1
     &           ,iphi0,iphi1
     &           ,n
     &           ,h
     &           ,ihlo0,ihlo1
     &           ,ihhi0,ihhi1
     &           ,m
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      REAL_T b(
     &           iblo0:ibhi0,
     &           iblo1:ibhi1)
      integer iplo0,iplo1
      integer iphi0,iphi1
      REAL_T p(
     &           iplo0:iphi0,
     &           iplo1:iphi1)
      REAL_T n
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL_T h(
     &           ihlo0:ihhi0,
     &           ihlo1:ihhi1)
      REAL_T m
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      Real_T pn
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      pn = min(one,p(i0,i1))**n
      a(i0,i1) = (pn * a(i0,i1)
     &     + (one - pn) *  b(i0,i1))
     &     * h(i0,i1)**m
      
      enddo
      enddo
      return
      end
      subroutine WATERDEPTH(
     &           depth
     &           ,idepthlo0,idepthlo1
     &           ,idepthhi0,idepthhi1
     &           ,thk
     &           ,ithklo0,ithklo1
     &           ,ithkhi0,ithkhi1
     &           ,usrf
     &           ,iusrflo0,iusrflo1
     &           ,iusrfhi0,iusrfhi1
     &           ,topg
     &           ,itopglo0,itopglo1
     &           ,itopghi0,itopghi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idepthlo0,idepthlo1
      integer idepthhi0,idepthhi1
      REAL_T depth(
     &           idepthlo0:idepthhi0,
     &           idepthlo1:idepthhi1)
      integer ithklo0,ithklo1
      integer ithkhi0,ithkhi1
      REAL_T thk(
     &           ithklo0:ithkhi0,
     &           ithklo1:ithkhi1)
      integer iusrflo0,iusrflo1
      integer iusrfhi0,iusrfhi1
      REAL_T usrf(
     &           iusrflo0:iusrfhi0,
     &           iusrflo1:iusrfhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL_T topg(
     &           itopglo0:itopghi0,
     &           itopglo1:itopghi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      depth(i0,i1) = (max(zero,
     &     usrf(i0,i1)
     &     -thk(i0,i1)
     &     -topg(i0,i1)))
      
      enddo
      enddo
      return
      end
      subroutine ZEROIFLESS(
     &           fab
     &           ,ifablo0,ifablo1
     &           ,ifabhi0,ifabhi1
     &           ,infab
     &           ,iinfablo0,iinfablo1
     &           ,iinfabhi0,iinfabhi1
     &           ,tol
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL_T fab(
     &           ifablo0:ifabhi0,
     &           ifablo1:ifabhi1)
      integer iinfablo0,iinfablo1
      integer iinfabhi0,iinfabhi1
      REAL_T infab(
     &           iinfablo0:iinfabhi0,
     &           iinfablo1:iinfabhi1)
      REAL_T tol
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      if (infab(i0,i1).lt.tol) then
         fab(i0,i1) = zero
      end if
      
      enddo
      enddo
      return
      end
      subroutine PWLFILL(
     &           a
     &           ,ialo0,ialo1
     &           ,iahi0,iahi1
     &           ,x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,xn
     &           ,ixnhi0
     &           ,bn
     &           ,ibnhi0
     &           ,dx
     &           ,idxhi0
     &           ,db
     &           ,idbhi0
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1)
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL_T x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1)
      integer ixnhi0
      REAL_T xn(
     &           0:ixnhi0)
      integer ibnhi0
      REAL_T bn(
     &           0:ibnhi0)
      integer idxhi0
      REAL_T dx(
     &           0:idxhi0)
      integer idbhi0
      REAL_T db(
     &           0:idbhi0)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, m, n
      Real_T xi
      do m = 0, ixnhi0 - 1
         dx(m) = xn(m+1) - xn(m)
         db(m) = bn(m+1) - bn(m)
      end do
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      xi = x(i0,i1)
      if (xi .lt. xn(0)) then
         a(i0,i1) = bn(0)
      else if (xi .ge. xn(ixnhi0)) then
         a(i0,i1) = bn(ixnhi0)
      else
         do m = 0, ixnhi0-1
            if ( xi.lt.xn(m+1) ) then
               a(i0,i1) = bn(m) + db(m)* (xi - xn(m))/dx(m)
               exit
            end if
         end do
      end if
      
      enddo
      enddo
      return
      end
      subroutine ABSLIMITFAB(
     &           u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,limit
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      REAL_T limit
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, n
      REAL_T abslim
      abslim = limit
      do n  = 0, nucomp-1
         
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         if (u(i0,i1,n).lt.zero) then
            u(i0,i1,n) = max(u(i0,i1,n),-abslim)
         else
            u(i0,i1,n) = min(u(i0,i1,n),abslim)
         end if
         
      enddo
      enddo
      end do
      return
      end
      subroutine MODLIMITFAB(
     &           u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,nlimit
     &           ,umax
     &           ,limit
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer nlimit
      REAL_T umax
      REAL_T limit
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, n
      REAL_T umod
      umax = zero
      nlimit = 0
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      umod = zero
      do n = 0, nucomp-1
         umod = umod + u(i0,i1,n)*u(i0,i1,n)
      end do
      umod = sqrt(max(zero,umod))
      umax = max(umax,umod)
      if (umod .gt. limit) then
         nlimit = nlimit + 1
         do n = 0, nucomp-1
            u(i0,i1,n) = limit * u(i0,i1,n)/umod
         end do
      end if
      
      enddo
      enddo
      return
      end
      subroutine MINFAB1(
     &           u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,limit
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      REAL_T limit
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1,n
      do n  = 0, nucomp-1
         
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         u(i0,i1,n) = min(u(i0,i1,n),limit)
         
      enddo
      enddo
      end do
      return
      end      
      subroutine MAXFAB1(
     &           u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,limit
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      REAL_T limit
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1,n
      do n  = 0, nucomp-1
         
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         u(i0,i1,n) = max(u(i0,i1,n),limit)
         
      enddo
      enddo
      end do
      return
      end 
      subroutine CHECKCOEF(
     &           ok
     &           ,ioklo0,ioklo1
     &           ,iokhi0,iokhi1
     &           ,dir
     &           ,mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,alpha
     &           ,ialphalo0,ialphalo1
     &           ,ialphahi0,ialphahi1
     &           ,mumin
     &           ,alphamin
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ioklo0,ioklo1
      integer iokhi0,iokhi1
      integer ok(
     &           ioklo0:iokhi0,
     &           ioklo1:iokhi1)
      integer dir
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL_T mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer ialphalo0,ialphalo1
      integer ialphahi0,ialphahi1
      REAL_T alpha(
     &           ialphalo0:ialphahi0,
     &           ialphalo1:ialphahi1)
      REAL_T mumin
      REAL_T alphamin
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ii0,ii1
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      if (alpha(i0,i1).gt.alphamin) then
         ok(i0,i1) = 1
      else if (mu(i0,i1).gt.mumin) then
         ok(i0,i1) = 1
      else if (mu(i0+ii0,i1+ii1).gt.mumin) then
         ok(i0,i1) = 1
      end if
      
      enddo
      enddo
      return
      end 
      subroutine APPLYMINVISCOSITY(
     &           mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,dir
     &           ,faceH
     &           ,ifaceHlo0,ifaceHlo1
     &           ,ifaceHhi0,ifaceHhi1
     &           ,ok
     &           ,ioklo0,ioklo1
     &           ,iokhi0,iokhi1
     &           ,mumin
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL_T mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer dir
      integer ifaceHlo0,ifaceHlo1
      integer ifaceHhi0,ifaceHhi1
      REAL_T faceH(
     &           ifaceHlo0:ifaceHhi0,
     &           ifaceHlo1:ifaceHhi1)
      integer ioklo0,ioklo1
      integer iokhi0,iokhi1
      integer ok(
     &           ioklo0:iokhi0,
     &           ioklo1:iokhi1)
      REAL_T mumin
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ii0,ii1
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
       
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      if (ok(i0,i1).ne.1) then
         if (faceH(i0,i1).gt.1.0d-10) then
            mu(i0,i1) = 
     &           max(mu(i0,i1), mumin);
         end if
         if (faceH(i0+ii0,i1+ii1).gt.1.0d-10) then
            mu(i0+ii0,i1+ii1) = 
     &           max(mu(i0+ii0,i1+ii1), mumin);
         end if
      end if
      
      enddo
      enddo
      return
      end
      subroutine UPWINDFLUXB(
     &           flux
     &           ,ifluxlo0,ifluxlo1
     &           ,ifluxhi0,ifluxhi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,s
     &           ,islo0,islo1
     &           ,ishi0,ishi1
     &           ,dir
     &           ,ifboxlo0,ifboxlo1
     &           ,ifboxhi0,ifboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifluxlo0,ifluxlo1
      integer ifluxhi0,ifluxhi1
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1)
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1)
      integer islo0,islo1
      integer ishi0,ishi1
      REAL_T s(
     &           islo0:ishi0,
     &           islo1:ishi1)
      integer dir
      integer ifboxlo0,ifboxlo1
      integer ifboxhi0,ifboxhi1
      integer i0,i1
      integer ii0,ii1
      REAL_T upws
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      
      do i1 = ifboxlo1,ifboxhi1
      do i0 = ifboxlo0,ifboxhi0

      if ( u(i0,i1) .gt. zero) then
         upws = s(i0-ii0,i1-ii1)
      else
         upws = s(i0,i1)
      end if
      flux(i0,i1) = upws * u(i0,i1)
      
      enddo
      enddo
      return
      end
      subroutine HOVERUNDERFLOTATION(
     &           hab
     &           ,ihablo0,ihablo1
     &           ,ihabhi0,ihabhi1
     &           ,hmin
     &           ,hmax
     &           ,thck
     &           ,ithcklo0,ithcklo1
     &           ,ithckhi0,ithckhi1
     &           ,topg
     &           ,itopglo0,itopglo1
     &           ,itopghi0,itopghi1
     &           ,rhoi
     &           ,rhoo
     &           ,sealevel
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ihablo0,ihablo1
      integer ihabhi0,ihabhi1
      REAL_T hab(
     &           ihablo0:ihabhi0,
     &           ihablo1:ihabhi1)
      REAL_T hmin
      REAL_T hmax
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL_T thck(
     &           ithcklo0:ithckhi0,
     &           ithcklo1:ithckhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL_T topg(
     &           itopglo0:itopghi0,
     &           itopglo1:itopghi1)
      REAL_T rhoi
      REAL_T rhoo
      REAL_T sealevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T ratio
      ratio = rhoo / rhoi
      hmax = -1.0e+12
      hmin = +1.0e+12
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      hab(i0,i1) = thck(i0,i1) - 
     &     ratio*max(zero, seaLevel - topg(i0,i1))
      hmax = max( hab(i0,i1), hmax)
      hmin = min( hab(i0,i1), hmin)
      
      enddo
      enddo
      return
      end
      subroutine INTEGRATEHEAVISIDE2D(
     &           r
     &           ,irlo0,irlo1
     &           ,irhi0,irhi1
     &           ,h
     &           ,ihlo0,ihlo1
     &           ,ihhi0,ihhi1
     &           ,n
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer irlo0,irlo1
      integer irhi0,irhi1
      REAL_T r(
     &           irlo0:irhi0,
     &           irlo1:irhi1)
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL_T h(
     &           ihlo0:ihhi0,
     &           ihlo1:ihhi1)
      integer n
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ix0,ix1
      integer iy0,iy1
      integer ixy0,ixy1
      integer kx,ky,nk,sx,sy
      REAL_T h00,h01,h10,h11,hh,rr,xx,yy,subdx,subdx2
      REAL_T w00(0:n-1,0:n-1),w01(0:n-1,0:n-1)
      REAL_T w10(0:n-1,0:n-1),w11(0:n-1,0:n-1)
      ix0 = CHF_ID(0,0)
                ix1 = CHF_ID(0,1)
      iy0 = CHF_ID(1,0)
                iy1 = CHF_ID(1,1)   
      ixy0 = CHF_ID(1,0)
                ixy1 = CHF_ID(1,1)   
#if (CH_SPACEDIM == 1)
      write(*,*) 'integrateheaviside2D called in 1D!'
      call MAYDAYERROR()
#elif (CH_SPACEDIM == 2)
      subdx = one / (n*two)
      subdx2 = subdx**2
      do ky = 0,n-1
         do kx = 0,n-1
           xx = (kx + half)*subdx
           yy = (ky + half)*subdx
           w00(kx,ky) = (one-xx)*(one-yy)
           w01(kx,ky) = (one-xx)*yy
           w10(kx,ky) = xx*(one-yy)
           w11(kx,ky) = xx*yy
        end do
      end do
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      rr = zero
      do sx = -1,1,2
         do sy = -1,1,2
            h00 = h(i0,i1)
            h10 = h(i0+sx,i1)
            h01 = h(i0,i1+sy)
            h11 = h(i0+sx,i1+sy)
            do ky = 0,n-1  
               do kx = 0,n-1
                  hh = h00*w01(kx,ky) 
     &                 + h01*w01(kx,ky) 
     &                 + h10*w10(kx,ky)
     &                 + h11*w11(kx,ky)
                  if (hh .gt. zero) then
                     rr = rr + subdx2
                  end if
               end do
            end do
         end do
      end do
      r(i0,i1) = rr
      
      enddo
      enddo
#endif
      return
      end
