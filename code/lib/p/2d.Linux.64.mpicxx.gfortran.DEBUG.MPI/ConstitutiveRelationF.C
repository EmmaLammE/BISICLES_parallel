#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine STRAININVARSSA(
     &           epsSqr
     &           ,iepsSqrlo0,iepsSqrlo1
     &           ,iepsSqrhi0,iepsSqrhi1
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
      integer iepsSqrlo0,iepsSqrlo1
      integer iepsSqrhi0,iepsSqrhi1
      REAL_T epsSqr(
     &           iepsSqrlo0:iepsSqrhi0,
     &           iepsSqrlo1:iepsSqrhi1)
      integer nderivscomp
      integer iderivslo0,iderivslo1
      integer iderivshi0,iderivshi1
      REAL_T derivs(
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
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

#if (CH_SPACEDIM == 1)
        epsSqr(i0,i1) = derivs(i0,i1,dudxcomp)**2
#elif (CH_SPACEDIM == 2)
        epsSqr(i0,i1) = derivs(i0,i1,dudxcomp)**2
     &                         +derivs(i0,i1,dvdycomp)**2
     &                         +(derivs(i0,i1,dudxcomp)
     &                          *derivs(i0,i1,dvdycomp))
     &                    +0.25*(derivs(i0,i1,dudycomp)
     &                          +derivs(i0,i1,dvdxcomp))**2
#endif 
      
      enddo
      enddo
      return 
      end
      subroutine COMPUTEARRHENIUSA(
     &           A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,theta
     &           ,ithetalo0,ithetalo1
     &           ,ithetahi0,ithetahi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,nParameter
     &           ,mParameter
     &           ,B0
     &           ,thetaR
     &           ,Kexponent
     &           ,Cfactor
     &           ,RgasConst
     &           ,Q
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL_T A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer ithetalo0,ithetalo1
      integer ithetahi0,ithetahi1
      REAL_T theta(
     &           ithetalo0:ithetahi0,
     &           ithetalo1:ithetahi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T nParameter
      REAL_T mParameter
      REAL_T B0
      REAL_T thetaR
      REAL_T Kexponent
      REAL_T Cfactor
      REAL_T RgasConst
      REAL_T Q
      integer i0,i1
      REAL_T APrefactor, QonR, thisA
      APrefactor = mParameter*((one/B0)**nParameter)      
      QonR = Q/RgasConst
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      thisA = exp( (three*Cfactor/((thetaR - theta(i0,i1))**Kexponent))
     &     - QonR/theta(i0,i1))
      A(i0,i1) = APrefactor*thisA
      
      enddo
      enddo
      end subroutine
      subroutine COMPUTEPATERSONA(
     &           A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,theta
     &           ,ithetalo0,ithetalo1
     &           ,ithetahi0,ithetahi1
     &           ,theta0
     &           ,itheta0lo0,itheta0lo1
     &           ,itheta0hi0,itheta0hi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,E
     &           ,A0
     &           ,R
     &           ,Qm
     &           ,Qp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL_T A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer ithetalo0,ithetalo1
      integer ithetahi0,ithetahi1
      REAL_T theta(
     &           ithetalo0:ithetahi0,
     &           ithetalo1:ithetahi1)
      integer itheta0lo0,itheta0lo1
      integer itheta0hi0,itheta0hi1
      REAL_T theta0(
     &           itheta0lo0:itheta0hi0,
     &           itheta0lo1:itheta0hi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T E
      REAL_T A0
      REAL_T R
      REAL_T Qm
      REAL_T Qp
      integer i0,i1
      REAL_T EA0, QmOnR, QpOnR, thisA, QonR, T, T0
      EA0 = E * A0     
      QmOnR = Qm/R
      QpOnR = Qp/R
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      T = theta(i0,i1)
      T0 = theta0(i0,i1)
      if (T.ge.T0) then
         QonR = QpOnR
      else
         QonR = QmOnR
      end if
      A(i0,i1) = EA0 * exp(-QonR *(one/T - one/T0))
      
      enddo
      enddo
      end subroutine
      subroutine COMPUTEZWINGERA(
     &           A
     &           ,iAlo0,iAlo1
     &           ,iAhi0,iAhi1
     &           ,theta
     &           ,ithetalo0,ithetalo1
     &           ,ithetahi0,ithetahi1
     &           ,theta0
     &           ,itheta0lo0,itheta0lo1
     &           ,itheta0hi0,itheta0hi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,E
     &           ,Ap
     &           ,R
     &           ,Qm
     &           ,Qp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iAlo0,iAlo1
      integer iAhi0,iAhi1
      REAL_T A(
     &           iAlo0:iAhi0,
     &           iAlo1:iAhi1)
      integer ithetalo0,ithetalo1
      integer ithetahi0,ithetahi1
      REAL_T theta(
     &           ithetalo0:ithetahi0,
     &           ithetalo1:ithetahi1)
      integer itheta0lo0,itheta0lo1
      integer itheta0hi0,itheta0hi1
      REAL_T theta0(
     &           itheta0lo0:itheta0hi0,
     &           itheta0lo1:itheta0hi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T E
      REAL_T Ap
      REAL_T R
      REAL_T Qm
      REAL_T Qp
      integer i0,i1
      REAL_T  A0,Am,QmOnR, QpOnR, thisA, QonR, T, T0
      QmOnR = Qm/R
      QpOnR = Qp/R
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      T = theta(i0,i1)
      T0 = theta0(i0,i1)
      Am =  Ap * exp(-QpOnR *(one/T0)) / exp(-QmOnR *(one/T0))
      if (T.ge.T0) then
         QonR = QpOnR
         A0 = Ap
      else
         QonR = QmOnR
         A0 = Am
      end if
      A(i0,i1) = E * A0 * exp(-QonR *(one/T))
      
      enddo
      enddo
      end subroutine
      subroutine COMPUTEGLENSMU0(
     &           mu0
     &           ,imu0lo0,imu0lo1
     &           ,imu0hi0,imu0hi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,nParameter
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imu0lo0,imu0lo1
      integer imu0hi0,imu0hi1
      REAL_T mu0(
     &           imu0lo0:imu0hi0,
     &           imu0lo1:imu0hi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T nParameter
      integer i0,i1
      REAL_T  muexponent
      muexponent = -1.0 / nParameter
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      mu0(i0,i1) =  half * (mu0(i0,i1)**muExponent)
      
      enddo
      enddo
      return
      end 
      subroutine COMPUTEGLENSMU(
     &           mu
     &           ,imulo0,imulo1
     &           ,imuhi0,imuhi1
     &           ,epsSqr
     &           ,iepsSqrlo0,iepsSqrlo1
     &           ,iepsSqrhi0,iepsSqrhi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,nExponent
     &           ,epsSqr0
     &           ,delta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imulo0,imulo1
      integer imuhi0,imuhi1
      REAL_T mu(
     &           imulo0:imuhi0,
     &           imulo1:imuhi1)
      integer iepsSqrlo0,iepsSqrlo1
      integer iepsSqrhi0,iepsSqrhi1
      REAL_T epsSqr(
     &           iepsSqrlo0:iepsSqrhi0,
     &           iepsSqrlo1:iepsSqrhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL_T nExponent
      REAL_T epsSqr0
      REAL_T delta
      integer i0,i1
      REAL_T epsExponent
      REAL_T thisMu, thisEpsSqr
      epsExponent = (one - nExponent)/(two*nExponent)
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

        thisEpsSqr = max(epsSqr(i0,i1),epsSqr0)
        thisMu=mu(i0,i1)*((thisEpsSqr)**epsExponent + delta)
        mu(i0,i1) = thisMu
      
      enddo
      enddo
      return
      end
