#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine UPLUSUT(
     &           e
     &           ,ielo0,ielo1
     &           ,iehi0,iehi1
     &           ,necomp
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,xxComp
     &           ,xyComp
     &           ,yxComp
     &           ,yyComp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer necomp
      integer ielo0,ielo1
      integer iehi0,iehi1
      REAL_T e(
     &           ielo0:iehi0,
     &           ielo1:iehi1,
     &           0:necomp-1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer xxComp
      integer xyComp
      integer yxComp
      integer yyComp
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T uxx,uxy,uyx,uyy
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if (CH_SPACEDIM == 1)
      uxx = u(i0,i1,xxComp)
      e(i0,i1,0) = uxx
#elif (CH_SPACEDIM == 2)
      uxx = u(i0,i1,xxComp)
      uyy = u(i0,i1,yyComp)
      uxy = u(i0,i1,xyComp)
      uyx = u(i0,i1,yxComp)
      e(i0,i1,xxComp) = uxx;
      e(i0,i1,xyComp) = half * ( uxy + uyx ) ;
      e(i0,i1,yxComp) = e(i0,i1,xyComp);
      e(i0,i1,yyComp) = uyy;
#endif
      
      enddo
      enddo
      return 
      end
      subroutine EPLUSTRE(
     &           e
     &           ,ielo0,ielo1
     &           ,iehi0,iehi1
     &           ,necomp
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,xxComp
     &           ,xyComp
     &           ,yxComp
     &           ,yyComp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer necomp
      integer ielo0,ielo1
      integer iehi0,iehi1
      REAL_T e(
     &           ielo0:iehi0,
     &           ielo1:iehi1,
     &           0:necomp-1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer xxComp
      integer xyComp
      integer yxComp
      integer yyComp
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T uxx,uxy,uyx,uyy,uzz
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if (CH_SPACEDIM == 1)
      e(i0,i1,0) = 2.0d0*u(i0,i1,0)
#elif (CH_SPACEDIM == 2)
      uxx = u(i0,i1,xxComp)
      uyy = u(i0,i1,yyComp)
      uxy = u(i0,i1,xyComp)
      uyx = u(i0,i1,yxComp)
      uzz = uxx + uyy
      e(i0,i1,xxComp) = (uxx + uzz) ;
      e(i0,i1,xyComp) = half * ( uxy + uyx ) ;
      e(i0,i1,yxComp) = e(i0,i1,xyComp);
      e(i0,i1,yyComp) = (uyy + uzz);
#endif
      
      enddo
      enddo
      return 
      end
      subroutine SYMTEIGEN(
     &           lambda
     &           ,ilambdalo0,ilambdalo1
     &           ,ilambdahi0,ilambdahi1
     &           ,nlambdacomp
     &           ,T
     &           ,iTlo0,iTlo1
     &           ,iThi0,iThi1
     &           ,nTcomp
     &           ,xxComp
     &           ,xyComp
     &           ,yxComp
     &           ,yyComp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           0:nlambdacomp-1)
      integer nTcomp
      integer iTlo0,iTlo1
      integer iThi0,iThi1
      REAL_T T(
     &           iTlo0:iThi0,
     &           iTlo1:iThi1,
     &           0:nTcomp-1)
      integer xxComp
      integer xyComp
      integer yxComp
      integer yyComp
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T b,d,txx,txy,tyy
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if (CH_SPACEDIM == 1)
      lambda(i0,i1,0) = T(i0,i1,0)
#elif (CH_SPACEDIM == 2)
      txx = T(i0,i1,xxComp)
      tyy = T(i0,i1,yyComp)
      txy = half*(abs(T(i0,i1,xyComp)) + abs(T(i0,i1,yxComp)))
      b = half * (txx + tyy)
      d = ( (half*(txx - tyy))**2 + txy**2)**half
      lambda(i0,i1,0) = b + d;
      lambda(i0,i1,1) = b - d;
#endif
      
      enddo
      enddo
      return 
      end
      subroutine NYECREVASSEDEPTHTP(
     &           dfab
     &           ,idfablo0,idfablo1
     &           ,idfabhi0,idfabhi1
     &           ,dwfab
     &           ,idwfablo0,idwfablo1
     &           ,idwfabhi0,idwfabhi1
     &           ,hfab
     &           ,ihfablo0,ihfablo1
     &           ,ihfabhi0,ihfabhi1
     &           ,tpfab
     &           ,itpfablo0,itpfablo1
     &           ,itpfabhi0,itpfabhi1
     &           ,habfab
     &           ,ihabfablo0,ihabfablo1
     &           ,ihabfabhi0,ihabfabhi1
     &           ,rhoi
     &           ,rhoo
     &           ,g
     &           ,eps
     &           ,dmax
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idfablo0,idfablo1
      integer idfabhi0,idfabhi1
      REAL_T dfab(
     &           idfablo0:idfabhi0,
     &           idfablo1:idfabhi1)
      integer idwfablo0,idwfablo1
      integer idwfabhi0,idwfabhi1
      REAL_T dwfab(
     &           idwfablo0:idwfabhi0,
     &           idwfablo1:idwfabhi1)
      integer ihfablo0,ihfablo1
      integer ihfabhi0,ihfabhi1
      REAL_T hfab(
     &           ihfablo0:ihfabhi0,
     &           ihfablo1:ihfabhi1)
      integer itpfablo0,itpfablo1
      integer itpfabhi0,itpfabhi1
      REAL_T tpfab(
     &           itpfablo0:itpfabhi0,
     &           itpfablo1:itpfabhi1)
      integer ihabfablo0,ihabfablo1
      integer ihabfabhi0,ihabfabhi1
      REAL_T habfab(
     &           ihabfablo0:ihabfabhi0,
     &           ihabfablo1:ihabfabhi1)
      REAL_T rhoi
      REAL_T rhoo
      REAL_T g
      REAL_T eps
      REAL_T dmax
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T h,d,tp,hab,Q,Qw,ds,dtot,dw
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      h = hfab(i0,i1)
      tp = tpfab( i0,i1 )
      hab = habfab(i0,i1)
      dw = dwfab(i0,i1)
      Q = rhoi * g * h * h
      Qw = rhoo * g * h * h
      dtot  = ( - (dw * rhoo + hab*rhoi)*Q
     &     + rhoo * (Qw * dw + h * tp) )
     &     / ( (rhoo-rhoi)*Q + tp*rhoo + eps )
      ds  = (Qw * dw + h * tp)/(eps + Q + tp) 
      dfab(i0,i1)  = max(zero,min(max(ds,dtot),h*dmax))
      
      enddo
      enddo
      return 
      end
      subroutine NYECREVASSEDEPTHT(
     &           dfab
     &           ,idfablo0,idfablo1
     &           ,idfabhi0,idfabhi1
     &           ,dwfab
     &           ,idwfablo0,idwfablo1
     &           ,idwfabhi0,idwfabhi1
     &           ,hfab
     &           ,ihfablo0,ihfablo1
     &           ,ihfabhi0,ihfabhi1
     &           ,tfab
     &           ,itfablo0,itfablo1
     &           ,itfabhi0,itfabhi1
     &           ,habfab
     &           ,ihabfablo0,ihabfablo1
     &           ,ihabfabhi0,ihabfabhi1
     &           ,rhoi
     &           ,rhoo
     &           ,g
     &           ,eps
     &           ,dlim
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idfablo0,idfablo1
      integer idfabhi0,idfabhi1
      REAL_T dfab(
     &           idfablo0:idfabhi0,
     &           idfablo1:idfabhi1)
      integer idwfablo0,idwfablo1
      integer idwfabhi0,idwfabhi1
      REAL_T dwfab(
     &           idwfablo0:idwfabhi0,
     &           idwfablo1:idwfabhi1)
      integer ihfablo0,ihfablo1
      integer ihfabhi0,ihfabhi1
      REAL_T hfab(
     &           ihfablo0:ihfabhi0,
     &           ihfablo1:ihfabhi1)
      integer itfablo0,itfablo1
      integer itfabhi0,itfabhi1
      REAL_T tfab(
     &           itfablo0:itfabhi0,
     &           itfablo1:itfabhi1)
      integer ihabfablo0,ihabfablo1
      integer ihabfabhi0,ihabfabhi1
      REAL_T habfab(
     &           ihabfablo0:ihabfabhi0,
     &           ihabfablo1:ihabfabhi1)
      REAL_T rhoi
      REAL_T rhoo
      REAL_T g
      REAL_T eps
      REAL_T dlim
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T h,d,t,hab,ds,db,dw
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      h = hfab(i0,i1)
      t = tfab( i0,i1 )/(eps + hfab(i0,i1))
      hab = habfab(i0,i1)
      dw = dwfab(i0,i1)
      db = rhoi/(rhoo-rhoi) * ( t/(rhoi * g) -hab)
      ds  = t/(rhoi * g) + rhoo/rhoi * dw
      dfab(i0,i1) =  max(zero,min(max(ds,ds+db),h*dlim))
      
      enddo
      enddo
      return 
      end
      subroutine FABMAX(
     &           a
     &           ,ialo0,ialo1
     &           ,iahi0,iahi1
     &           ,b
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
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
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      a(i0,i1) = max( a(i0,i1), b(i0,i1))
      
      enddo
      enddo
      return 
      end
      subroutine FABMIN(
     &           a
     &           ,ialo0,ialo1
     &           ,iahi0,iahi1
     &           ,b
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
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
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      a(i0,i1) = min( a(i0,i1), b(i0,i1))
      
      enddo
      enddo
      return 
      end
      subroutine CREVASSEMU(
     &           mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,thck
     &           ,ithcklo0,ithcklo1
     &           ,ithckhi0,ithckhi1
     &           ,depth
     &           ,idepthlo0,idepthlo1
     &           ,idepthhi0,idepthhi1
     &           ,eps
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL_T mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL_T thck(
     &           ithcklo0:ithckhi0,
     &           ithcklo1:ithckhi1)
      integer idepthlo0,idepthlo1
      integer idepthhi0,idepthhi1
      REAL_T depth(
     &           idepthlo0:idepthhi0,
     &           idepthlo1:idepthhi1)
      REAL_T eps
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T h,d
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      d = depth(i0,i1)
      h = thck(i0,i1)
      mu(i0,i1) = mu(i0,i1)
     &    * h * (h - min(d,h-eps))/(h*h + eps)
      
      enddo
      enddo
      return 
      end
