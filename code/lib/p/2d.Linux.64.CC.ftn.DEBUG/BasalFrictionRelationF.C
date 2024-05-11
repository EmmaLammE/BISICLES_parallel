#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
#include "IceConstants.H"
      subroutine BFRICTIONPOWER(
     &           alpha
     &           ,ialphalo0,ialphalo1
     &           ,ialphahi0,ialphahi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,C
     &           ,iClo0,iClo1
     &           ,iChi0,iChi1
     &           ,p
     &           ,iplo0,iplo1
     &           ,iphi0,iphi1
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,m
     &           ,n
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialphalo0,ialphalo1
      integer ialphahi0,ialphahi1
      REAL_T alpha(
     &           ialphalo0:ialphahi0,
     &           ialphalo1:ialphahi1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer iClo0,iClo1
      integer iChi0,iChi1
      REAL_T C(
     &           iClo0:iChi0,
     &           iClo1:iChi1)
      integer iplo0,iplo1
      integer iphi0,iphi1
      REAL_T p(
     &           iplo0:iphi0,
     &           iplo1:iphi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1)
      REAL_T m
      REAL_T n
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T usq, mExp, nExp, usq0
      usq0 = TINY_VEL**2
      mExp = m/two
      nExp = n
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         
         usq = u(i0,i1,0)**2
     &        + u(i0,i1,1)**2
         if (mask(i0,i1).eq.GROUNDEDMASKVAL) then       
            alpha(i0,i1) = C(i0,i1) 
     &           * (usq + usq0)**mExp
            if (n.ne.0.0d0) then
               alpha(i0,i1) = alpha(i0,i1) *
     &            p(i0,i1)**nExp 
            end if
         else
            alpha(i0,i1) = C(i0,i1) 
         end if
      
      enddo
      enddo
      return
      end
      subroutine BFRICTIONAUU(
     &           alpha
     &           ,ialphalo0,ialphalo1
     &           ,ialphahi0,ialphahi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialphalo0,ialphalo1
      integer ialphahi0,ialphahi1
      REAL_T alpha(
     &           ialphalo0:ialphahi0,
     &           ialphalo1:ialphahi1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, ncomp, icomp
      REAL_T usq
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      
      usq = u(i0,i1,0)**2
     &     + u(i0,i1,1)**2
      alpha(i0,i1) = alpha(i0,i1) * usq
      
      enddo
      enddo
      return
      end
      subroutine BFRICTIONPLIMITTSAI(
     &           alpha
     &           ,ialphalo0,ialphalo1
     &           ,ialphahi0,ialphahi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,p
     &           ,iplo0,iplo1
     &           ,iphi0,iphi1
     &           ,a
     &           ,ialo0,ialo1
     &           ,iahi0,iahi1
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialphalo0,ialphalo1
      integer ialphahi0,ialphahi1
      REAL_T alpha(
     &           ialphalo0:ialphahi0,
     &           ialphalo1:ialphahi1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer iplo0,iplo1
      integer iphi0,iphi1
      REAL_T p(
     &           iplo0:iphi0,
     &           iplo1:iphi1)
      integer ialo0,ialo1
      integer iahi0,iahi1
      REAL_T a(
     &           ialo0:iahi0,
     &           ialo1:iahi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T usq, usq0
      usq0 = TINY_VEL**2
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      
      usq = usq0 + u(i0,i1,0)**2
     &     + u(i0,i1,1)**2
      alpha(i0,i1) = min(alpha(i0,i1),
     &     a(i0,i1)  * p(i0,i1) * usq**(-half))
      
      enddo
      enddo
      return
      end
      subroutine BFRICTIONPLIMITLEGUY(
     &           alpha
     &           ,ialphalo0,ialphalo1
     &           ,ialphahi0,ialphahi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,p
     &           ,iplo0,iplo1
     &           ,iphi0,iphi1
     &           ,a
     &           ,n
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialphalo0,ialphalo1
      integer ialphahi0,ialphahi1
      REAL_T alpha(
     &           ialphalo0:ialphahi0,
     &           ialphalo1:ialphahi1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      integer iplo0,iplo1
      integer iphi0,iphi1
      REAL_T p(
     &           iplo0:iphi0,
     &           iplo1:iphi1)
      REAL_T a
      REAL_T n
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T usq, usq0, pn, oneOnN
      usq0 = TINY_VEL**2
      oneOnN = one / n
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      
      usq = usq0 + u(i0,i1,0)**2
     &     + u(i0,i1,1)**2
      pn = p(i0,i1)**n
      alpha(i0,i1) = alpha(i0,i1)
     &     * ((pn) / (a*usq**(half) + pn))**oneOnN
      
      enddo
      enddo
      return
      end
      subroutine BFRICTIONJOUGHIN(
     &           alpha
     &           ,ialphalo0,ialphalo1
     &           ,ialphahi0,ialphahi1
     &           ,u
     &           ,iulo0,iulo1
     &           ,iuhi0,iuhi1
     &           ,nucomp
     &           ,ur
     &           ,N
     &           ,iNlo0,iNlo1
     &           ,iNhi0,iNhi1
     &           ,Nr
     &           ,m
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ialphalo0,ialphalo1
      integer ialphahi0,ialphahi1
      REAL_T alpha(
     &           ialphalo0:ialphahi0,
     &           ialphalo1:ialphahi1)
      integer nucomp
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           0:nucomp-1)
      REAL_T ur
      integer iNlo0,iNlo1
      integer iNhi0,iNhi1
      REAL_T N(
     &           iNlo0:iNhi0,
     &           iNlo1:iNhi1)
      REAL_T Nr
      REAL_T m
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T usq, usq0, q, om
      om = one / m
      usq0 = TINY_VEL**2
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      
      usq = usq0 + u(i0,i1,0)**2
     &     + u(i0,i1,1)**2
      q = (usq**half)/(abs(ur) + TINY_VEL)
      if (Nr .lt. TINY_THICKNESS) then
         alpha(i0,i1) = alpha(i0,i1)
     &        * ( q + one )**(-m)
      else
         alpha(i0,i1) = alpha(i0,i1)
     &        *     N(i0,i1)
     &        * ( q * Nr**om + N(i0,i1)**om )**(-m)
      end if
      
      enddo
      enddo
      return
      end
      subroutine BFRICTIONLEGUYEFFPRES(
     &           N
     &           ,iNlo0,iNlo1
     &           ,iNhi0,iNhi1
     &           ,hab
     &           ,ihablo0,ihablo1
     &           ,ihabhi0,ihabhi1
     &           ,h
     &           ,ihlo0,ihlo1
     &           ,ihhi0,ihhi1
     &           ,p
     &           ,rhog
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iNlo0,iNlo1
      integer iNhi0,iNhi1
      REAL_T N(
     &           iNlo0:iNhi0,
     &           iNlo1:iNhi1)
      integer ihablo0,ihablo1
      integer ihabhi0,ihabhi1
      REAL_T hab(
     &           ihablo0:ihabhi0,
     &           ihablo1:ihabhi1)
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL_T h(
     &           ihlo0:ihhi0,
     &           ihlo1:ihhi1)
      REAL_T p
      REAL_T rhog
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL_T pp
      pp = one - p
      
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

      N(i0,i1) = rhog 
     &     * h(i0,i1)**pp
     &     * hab(i0,i1)**p
      
      enddo
      enddo
      return
      end
