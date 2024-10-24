      subroutine L1L2LIMITFAB(
     &           fab
     &           ,ifablo0,ifablo1
     &           ,ifabhi0,ifabhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,mv
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     &           ifablo0:ifabhi0,
     &           ifablo1:ifabhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 mv
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (.not.(fab(i0,i1) .lt. mv)) then
         fab(i0,i1) = mv
      end if
      enddo
      enddo
      return
      end subroutine
      subroutine L1L2SUPRESSFAB(
     &           fab
     &           ,ifablo0,ifablo1
     &           ,ifabhi0,ifabhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,mv
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     &           ifablo0:ifabhi0,
     &           ifablo1:ifabhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 mv
      integer i0,i1
      REAL*8 t
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      t = fab(i0,i1)
      if (t .gt. (0.0d0)) then
         if (.not.(fab(i0,i1) .le. mv)) then
            fab(i0,i1) = mv + log(t - mv + 1.0)
         end if
      else
         if (.not.(fab(i0,i1) .ge. -mv)) then
            fab(i0,i1) = -mv - log(-t - mv + 1.0)
         end if 
      end if
      enddo
      enddo
      return
      end subroutine
      subroutine ANALYTICL1L2MU(
     &           mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,epsSqr
     &           ,iepsSqrlo0,iepsSqrlo1
     &           ,iepsSqrhi0,iepsSqrhi1
     &           ,torSqr
     &           ,itorSqrlo0,itorSqrlo1
     &           ,itorSqrhi0,itorSqrhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,sigma
     &           ,nExponent
     &           ,epsSqr0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL*8 mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL*8 A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer iepsSqrlo0,iepsSqrlo1
      integer iepsSqrhi0,iepsSqrhi1
      REAL*8 epsSqr(
     &           iepsSqrlo0:iepsSqrhi0,
     &           iepsSqrlo1:iepsSqrhi1)
      integer itorSqrlo0,itorSqrlo1
      integer itorSqrhi0,itorSqrhi1
      REAL*8 torSqr(
     &           itorSqrlo0:itorSqrhi0,
     &           itorSqrlo1:itorSqrhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 sigma
      REAL*8 nExponent
      REAL*8 epsSqr0
      integer i0,i1, niter
      REAL*8  tmu, tA, tesqr, ttsqr, q, sigmaSqr  
      sigmaSqr = sigma*sigma;
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      tA = A(i0,i1)
      tesqr = epsSqr(i0,i1) + epsSqr0 
      ttsqr = sigmaSqr * torSqr(i0,i1)
      q = (4.0 * ta**(2.0d0) * ttsqr**(3.0d0) + 27.0 * tesqr)**(0.500d0)
     &     + (3.0d0)**((3.0d0)/(2.0d0)) * tesqr**(0.500d0) 
      q = q**((1.000d0 / 3.000d0))
      tmu = q**(2.0d0) - (2.0d0)**((2.0d0)/(3.0d0)) * tA**((2.0d0)/(3.0d
     &0)) * ttsqr
      tmu = tmu / ((4.0d0)**((2.0d0)/(3.0d0)) * (3.0d0)**(0.500d0) * tA*
     &*((1.000d0 / 3.000d0)) 
     &     * tesqr**(0.500d0) * q)
      mu(i0,i1) = tmu
      enddo
      enddo
      return 
      end
      subroutine COMPUTEL1L2RES(
     &           res
     &           ,mu
     &           ,A
     &           ,esqr
     &           ,tsqr
     &           ,p
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 res
      REAL*8 mu
      REAL*8 A
      REAL*8 esqr
      REAL*8 tsqr
      REAL*8 p
      res = 1.0 - 2.0 * A * mu * (4.0 * mu**2 * esqr + tsqr)**p
      return
      end subroutine
      subroutine L1L2SECANTSOLVE(
     &           mu
     &           ,mua
     &           ,A
     &           ,esqr
     &           ,tsqr
     &           ,p
     &           ,niter
     &           ,tol
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 mu
      REAL*8 mua
      REAL*8 A
      REAL*8 esqr
      REAL*8 tsqr
      REAL*8 p
      integer niter
      REAL*8 tol
      REAL*8 fa,fb,t
      integer i, j
      if (abs(mua - mu) .lt. 0.005) then
         mua = 1.005 * mu;
      end if
      call COMPUTEL1L2RES(fa,mua,A,esqr,tsqr,p)
      do i = 0, niter
         call COMPUTEL1L2RES(fb,mu,A,esqr,tsqr,p) 
         if (abs(fb) .lt. tol) exit
         t = mu
         mu = mua - (mua - mu)/(fa - fb)*fa
         mua = t
         fa = fb 
      end do
      mua = fb
      return
      end subroutine
      subroutine COMPUTEL1L2MU(
     &           mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,res
     &           ,ireslo0,ireslo1
     &           ,ireshi0,ireshi1
     &           ,A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,epsSqr
     &           ,iepsSqrlo0,iepsSqrlo1
     &           ,iepsSqrhi0,iepsSqrhi1
     &           ,torSqr
     &           ,itorSqrlo0,itorSqrlo1
     &           ,itorSqrhi0,itorSqrhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,sigma
     &           ,nExponent
     &           ,epsSqr0
     &           ,delta
     &           ,tol
     &           ,maxIter
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL*8 mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL*8 res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1)
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL*8 A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer iepsSqrlo0,iepsSqrlo1
      integer iepsSqrhi0,iepsSqrhi1
      REAL*8 epsSqr(
     &           iepsSqrlo0:iepsSqrhi0,
     &           iepsSqrlo1:iepsSqrhi1)
      integer itorSqrlo0,itorSqrlo1
      integer itorSqrhi0,itorSqrhi1
      REAL*8 torSqr(
     &           itorSqrlo0:itorSqrhi0,
     &           itorSqrlo1:itorSqrhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 sigma
      REAL*8 nExponent
      REAL*8 epsSqr0
      REAL*8 delta
      REAL*8 tol
      integer maxIter
      integer i0,i1, niter, isec
      REAL*8 tmua, tmu, tA, tesqr, ttsqr, p1, p2, p3 , mu0, t
      REAL*8 fa,fb
      p1 = (nExponent - (1.0d0))*(0.500d0)
      p2 = -p1 / nExponent
      p3 = - (1.0d0) / nExponent
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
          tA = A(i0,i1)
          tesqr = epsSqr(i0,i1) + epsSqr0 
          ttsqr = sigma * sigma * torSqr(i0,i1)
          if (ttsqr .gt. epsSqr0) then
             niter = maxIter
             tmu = mu(i0,i1)
             tmua = res(i0,i1)
             if (abs(tmua - tmu) .lt. 0.005d0*(tmu+1.0)) then
                tmua = 1.005d0 * tmu;
             end if       
             fa = 1.0d0 - 2.0d0 * tA * tmua * 
     &            (4.0d0 * tmua**2 * tesqr + ttsqr)**p1
             do isec = 0, niter
                fb = 1.0d0 - 2.0d0 * tA * tmu * 
     &               (4.0d0 * tmu**2 * tesqr + ttsqr)**p1
                if (abs(fb) .lt. tol) exit
                if (abs(fb-fa) .lt. tol**2) exit
                t = tmu
                tmu = tmua - (tmua - tmu)/(fa - fb)*fa 
                tmua = t
                fa = fb 
             end do
             tmua = fb
             res(i0,i1) = fb
             mu(i0,i1) = tmu
          else
             tmua = (0.500d0) * tA**p3 * tesqr**p2
             mu(i0,i1) = tmua
             res(i0,i1) = 0.0d0
          end if
          mu(i0,i1) = mu(i0,i1) + tA**p3 * delta
      enddo
      enddo
      return
      end subroutine
      subroutine L1L2PHIINTEGRANDS(
     &           sixx
     &           ,isixxlo0,isixxlo1
     &           ,isixxhi0,isixxhi1
     &           ,sixy
     &           ,isixylo0,isixylo1
     &           ,isixyhi0,isixyhi1
     &           ,siyy
     &           ,isiyylo0,isiyylo1
     &           ,isiyyhi0,isiyyhi1
     &           ,H
     &           ,iHlo0,iHlo1
     &           ,iHhi0,iHhi1
     &           ,mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,derivs
     &           ,iderivslo0,iderivslo1
     &           ,iderivshi0,iderivshi1
     &           ,nderivscomp
     &           ,dudxComp
     &           ,dudyComp
     &           ,dvdxComp
     &           ,dvdyComp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isixxlo0,isixxlo1
      integer isixxhi0,isixxhi1
      REAL*8 sixx(
     &           isixxlo0:isixxhi0,
     &           isixxlo1:isixxhi1)
      integer isixylo0,isixylo1
      integer isixyhi0,isixyhi1
      REAL*8 sixy(
     &           isixylo0:isixyhi0,
     &           isixylo1:isixyhi1)
      integer isiyylo0,isiyylo1
      integer isiyyhi0,isiyyhi1
      REAL*8 siyy(
     &           isiyylo0:isiyyhi0,
     &           isiyylo1:isiyyhi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     &           iHlo0:iHhi0,
     &           iHlo1:iHhi1)
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL*8 mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer nderivscomp
      integer iderivslo0,iderivslo1
      integer iderivshi0,iderivshi1
      REAL*8 derivs(
     &           iderivslo0:iderivshi0,
     &           iderivslo1:iderivshi1,
     &           0:nderivscomp-1)
      integer dudxComp
      integer dudyComp
      integer dvdxComp
      integer dvdyComp
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 thisMuH
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
             thisMuH = mu(i0,i1) * H(i0,i1)
             sixx(i0,i1) = thisMuH * 
     &            ((4.0d0) *  derivs(i0,i1,dudxcomp)
     &            + (2.0d0) *  derivs(i0,i1,dvdycomp))
             sixy(i0,i1) = thisMuH * 
     &            (derivs(i0,i1,dudycomp)
     &            + derivs(i0,i1,dvdxcomp))
             siyy(i0,i1) = thisMuH * 
     &            ((4.0d0) *  derivs(i0,i1,dvdycomp)
     &            + (2.0d0) *  derivs(i0,i1,dudxcomp))        
      enddo
      enddo
      return
      end subroutine
      subroutine L1L2UINTEGRANDS(
     &           fzx
     &           ,ifzxlo0,ifzxlo1
     &           ,ifzxhi0,ifzxhi1
     &           ,fzy
     &           ,ifzylo0,ifzylo1
     &           ,ifzyhi0,ifzyhi1
     &           ,H
     &           ,iHlo0,iHlo1
     &           ,iHhi0,iHhi1
     &           ,A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,epssqr
     &           ,iepssqrlo0,iepssqrlo1
     &           ,iepssqrhi0,iepssqrhi1
     &           ,phizx
     &           ,iphizxlo0,iphizxlo1
     &           ,iphizxhi0,iphizxhi1
     &           ,phizy
     &           ,iphizylo0,iphizylo1
     &           ,iphizyhi0,iphizyhi1
     &           ,n
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifzxlo0,ifzxlo1
      integer ifzxhi0,ifzxhi1
      REAL*8 fzx(
     &           ifzxlo0:ifzxhi0,
     &           ifzxlo1:ifzxhi1)
      integer ifzylo0,ifzylo1
      integer ifzyhi0,ifzyhi1
      REAL*8 fzy(
     &           ifzylo0:ifzyhi0,
     &           ifzylo1:ifzyhi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     &           iHlo0:iHhi0,
     &           iHlo1:iHhi1)
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL*8 A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL*8 mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer iepssqrlo0,iepssqrlo1
      integer iepssqrhi0,iepssqrhi1
      REAL*8 epssqr(
     &           iepssqrlo0:iepssqrhi0,
     &           iepssqrlo1:iepssqrhi1)
      integer iphizxlo0,iphizxlo1
      integer iphizxhi0,iphizxhi1
      REAL*8 phizx(
     &           iphizxlo0:iphizxhi0,
     &           iphizxlo1:iphizxhi1)
      integer iphizylo0,iphizylo1
      integer iphizyhi0,iphizyhi1
      REAL*8 phizy(
     &           iphizylo0:iphizyhi0,
     &           iphizylo1:iphizyhi1)
      REAL*8 n
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 phisqr, g , p
      p = (0.500d0) * (n-(1.0d0))
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
        phisqr = phizx(i0,i1)**2 
     &     +  phizy(i0,i1)**2
        g = (2.0d0) * A(i0,i1) * H(i0,i1)
     &       * ((4.0d0) * mu(i0,i1)**2 
     &       * epsSqr(i0,i1) + phisqr)**p
       fzx(i0,i1) = g * phizx(i0,i1)
       fzy(i0,i1) = g * phizy(i0,i1)
      enddo
      enddo
      return
      end subroutine
      subroutine L1L2UIGRAND(
     &           fzi
     &           ,ifzilo0,ifzilo1
     &           ,ifzihi0,ifzihi1
     &           ,nfzicomp
     &           ,H
     &           ,iHlo0,iHlo1
     &           ,iHhi0,iHhi1
     &           ,A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,epssqr
     &           ,iepssqrlo0,iepssqrlo1
     &           ,iepssqrhi0,iepssqrhi1
     &           ,phizi
     &           ,iphizilo0,iphizilo1
     &           ,iphizihi0,iphizihi1
     &           ,nphizicomp
     &           ,n
     &           ,ncomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfzicomp
      integer ifzilo0,ifzilo1
      integer ifzihi0,ifzihi1
      REAL*8 fzi(
     &           ifzilo0:ifzihi0,
     &           ifzilo1:ifzihi1,
     &           0:nfzicomp-1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     &           iHlo0:iHhi0,
     &           iHlo1:iHhi1)
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL*8 A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL*8 mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer iepssqrlo0,iepssqrlo1
      integer iepssqrhi0,iepssqrhi1
      REAL*8 epssqr(
     &           iepssqrlo0:iepssqrhi0,
     &           iepssqrlo1:iepssqrhi1)
      integer nphizicomp
      integer iphizilo0,iphizilo1
      integer iphizihi0,iphizihi1
      REAL*8 phizi(
     &           iphizilo0:iphizihi0,
     &           iphizilo1:iphizihi1,
     &           0:nphizicomp-1)
      REAL*8 n
      integer ncomp
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, icomp
      REAL*8 phisqr, g , p
      p = (0.500d0) * (n-(1.0d0))
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      phisqr = 0.0;
      do icomp = 0, ncomp -1
         phisqr = phisqr +  phizi(i0,i1,icomp)**2 
      end do
      g = (2.0d0) * min( A(i0,i1), 1.0e-10) * H(i0,i1)
     &     * ((4.0d0) * mu(i0,i1)**2 
     &     * epsSqr(i0,i1) + phisqr)**p
       do icomp = 0, ncomp -1
          fzi(i0,i1,icomp) = g * phizi(i0,i1,icomp)
       end do
      enddo
      enddo
      return
      end subroutine
      subroutine L1L2MODLIMIT(
     &           fab
     &           ,ifablo0,ifablo1
     &           ,ifabhi0,ifabhi1
     &           ,nfabcomp
     &           ,lim
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,ncomp
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfabcomp
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     &           ifablo0:ifabhi0,
     &           ifablo1:ifabhi1,
     &           0:nfabcomp-1)
      REAL*8 lim
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ncomp
      REAL*8 limsq,t
      integer i0,i1, icomp
      limsq = lim**(2.0d0);
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      t = (0.0d0)
      do icomp = 0, ncomp-1
         t = t + fab(i0,i1,icomp)**(2.0d0)
      end do
      if (t .gt. limsq) then
         t = lim * t**(-(0.500d0))
         do icomp = 0, ncomp-1
            fab(i0,i1,icomp) = t * fab(i0,i1,icomp) 
         end do
      end if
      enddo
      enddo
      return
      end
      subroutine L1L2ADDGSTRESS(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,G
     &           ,iGlo0,iGlo1
     &           ,iGhi0,iGhi1
     &           ,nGcomp
     &           ,H
     &           ,iHlo0,iHlo1
     &           ,iHhi0,iHhi1
     &           ,sigma
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,ncomp
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nGcomp
      integer iGlo0,iGlo1
      integer iGhi0,iGhi1
      REAL*8 G(
     &           iGlo0:iGhi0,
     &           iGlo1:iGhi1,
     &           0:nGcomp-1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     &           iHlo0:iHhi0,
     &           iHlo1:iHhi1)
      REAL*8 sigma
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ncomp
      integer i0,i1, icomp
      do icomp = 0, ncomp-1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         phi(i0,i1,icomp) = phi(i0,i1,icomp)
     &        - sigma * G (i0,i1,icomp) * H(i0,i1)
      enddo
      enddo
      end do
      return
      end subroutine
      subroutine L1L2PHITILDESQR(
     &           ptsqr
     &           ,iptsqrlo0,iptsqrlo1
     &           ,iptsqrhi0,iptsqrhi1
     &           ,H
     &           ,iHlo0,iHlo1
     &           ,iHhi0,iHhi1
     &           ,grads
     &           ,igradslo0,igradslo1
     &           ,igradshi0,igradshi1
     &           ,ngradscomp
     &           ,rhog
     &           ,ncomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iptsqrlo0,iptsqrlo1
      integer iptsqrhi0,iptsqrhi1
      REAL*8 ptsqr(
     &           iptsqrlo0:iptsqrhi0,
     &           iptsqrlo1:iptsqrhi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     &           iHlo0:iHhi0,
     &           iHlo1:iHhi1)
      integer ngradscomp
      integer igradslo0,igradslo1
      integer igradshi0,igradshi1
      REAL*8 grads(
     &           igradslo0:igradshi0,
     &           igradslo1:igradshi1,
     &           0:ngradscomp-1)
      REAL*8 rhog
      integer ncomp
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, icomp
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      ptsqr(i0,i1) =  (rhog * H(i0,i1) 
     &     * grads( i0,i1, 0))**2
      enddo
      enddo
      if (ncomp .gt. 1) then
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         ptsqr(i0,i1) = ptsqr(i0,i1) 
     &        + (rhog * H(i0,i1) 
     &        * grads( i0,i1, 1))**2
      enddo
      enddo
      end if
      return 
      end subroutine
      subroutine L1L2DIFFFACTOR(
     &           D
     &           ,iDlo0,iDlo1
     &           ,iDhi0,iDhi1
     &           ,H
     &           ,iHlo0,iHlo1
     &           ,iHhi0,iHhi1
     &           ,grads
     &           ,igradslo0,igradslo1
     &           ,igradshi0,igradshi1
     &           ,ngradscomp
     &           ,n
     &           ,rhog
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iDlo0,iDlo1
      integer iDhi0,iDhi1
      REAL*8 D(
     &           iDlo0:iDhi0,
     &           iDlo1:iDhi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     &           iHlo0:iHhi0,
     &           iHlo1:iHhi1)
      integer ngradscomp
      integer igradslo0,igradslo1
      integer igradshi0,igradshi1
      REAL*8 grads(
     &           igradslo0:igradshi0,
     &           igradslo1:igradshi1,
     &           0:ngradscomp-1)
      REAL*8 n
      REAL*8 rhog
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, ncomp
      REAL*8  G, f, p1, p2
      f = 2.0*rhog**n
      p1 = (n+(2.0d0))
      p2 = ((n-(1.0d0))/(2.0d0))
      ncomp = ngradscomp
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      G = grads(i0,i1,0)**2
      if (ncomp .gt. 1) then 
         G = G + grads(i0,i1,1)**2
      end if
      D(i0,i1) = f * H(i0,i1)**p1 * G**p2
      enddo
      enddo
      return
      end
      subroutine L1L2COMPUTEDIFFUSIVITY(
     &           D
     &           ,iDlo0,iDlo1
     &           ,iDhi0,iDhi1
     &           ,ux
     &           ,iuxlo0,iuxlo1
     &           ,iuxhi0,iuxhi1
     &           ,uy
     &           ,iuylo0,iuylo1
     &           ,iuyhi0,iuyhi1
     &           ,thck
     &           ,ithcklo0,ithcklo1
     &           ,ithckhi0,ithckhi1
     &           ,s
     &           ,islo0,islo1
     &           ,ishi0,ishi1
     &           ,dir
     &           ,dx
     &           ,dy
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iDlo0,iDlo1
      integer iDhi0,iDhi1
      REAL*8 D(
     &           iDlo0:iDhi0,
     &           iDlo1:iDhi1)
      integer iuxlo0,iuxlo1
      integer iuxhi0,iuxhi1
      REAL*8 ux(
     &           iuxlo0:iuxhi0,
     &           iuxlo1:iuxhi1)
      integer iuylo0,iuylo1
      integer iuyhi0,iuyhi1
      REAL*8 uy(
     &           iuylo0:iuyhi0,
     &           iuylo1:iuyhi1)
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     &           ithcklo0:ithckhi0,
     &           ithcklo1:ithckhi1)
      integer islo0,islo1
      integer ishi0,ishi1
      REAL*8 s(
     &           islo0:ishi0,
     &           islo1:ishi1)
      integer dir
      REAL*8 dx
      REAL*8 dy
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ix0,ix1
      integer iy0,iy1
      REAL*8 oneOnDx, uxc,uyc, uc, dsx, dsy, ds
      oneOnDx = (1.0d0)/dx
      ix0 = CHF_ID(0,0)
      ix1 = CHF_ID(0,1)
      iy0 = CHF_ID(1,0)
      iy1 = CHF_ID(1,1)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      uxc = (0.500d0) * (ux(i0,i1)+ux(i0+ix0,i1+ix1))
      uyc = (0.500d0) * (uy(i0,i1)+uy(i0+iy0,i1+iy1))
      uc = (uxc**2 + uyc**2)**(0.500d0)
      dsx = (0.500d0)*(s(i0+ix0,i1+ix1)
     &     -s(i0-ix0,i1-ix1))/dx
      dsy = (0.500d0)*(s(i0+iy0,i1+iy1)
     &     -s(i0-iy0,i1-iy1))/dy
      ds = (dsx**2 + dsy**2)**(0.500d0)
      D(i0,i1) = thck(i0,i1) * uc
     &		       / max(ds,1.0e-4)
      enddo
      enddo
      return 
      end
      subroutine L1L2CELLTOFACEHARMONIC(
     &           face
     &           ,ifacelo0,ifacelo1
     &           ,ifacehi0,ifacehi1
     &           ,cell
     &           ,icelllo0,icelllo1
     &           ,icellhi0,icellhi1
     &           ,dir
     &           ,tiny
     &           ,ifboxlo0,ifboxlo1
     &           ,ifboxhi0,ifboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacelo0,ifacelo1
      integer ifacehi0,ifacehi1
      REAL*8 face(
     &           ifacelo0:ifacehi0,
     &           ifacelo1:ifacehi1)
      integer icelllo0,icelllo1
      integer icellhi0,icellhi1
      REAL*8 cell(
     &           icelllo0:icellhi0,
     &           icelllo1:icellhi1)
      integer dir
      REAL*8 tiny
      integer ifboxlo0,ifboxlo1
      integer ifboxhi0,ifboxhi1
      integer i0,i1
      integer ii0,ii1
      ii0 = CHF_ID(dir,0)
      ii1 = CHF_ID(dir,1)
      do i1 = ifboxlo1,ifboxhi1
      do i0 = ifboxlo0,ifboxhi0
      face(i0,i1) = (0.500d0) * 
     &     cell(i0,i1) * cell (i0-ii0,i1-ii1) / 
     &     (tiny + cell(i0,i1) + cell (i0-ii0,i1-ii1))
      enddo
      enddo
      return
      end
