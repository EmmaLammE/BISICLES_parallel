#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine DIVERGENCE(
     &           uEdge
     &           ,iuEdgelo0,iuEdgelo1
     &           ,iuEdgehi0,iuEdgehi1
     &           ,nuEdgecomp
     &           ,div
     &           ,idivlo0,idivlo1
     &           ,idivhi0,idivhi1
     &           ,ndivcomp
     &           ,igridIntlo0,igridIntlo1
     &           ,igridInthi0,igridInthi1
     &           ,dx
     &           ,idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nuEdgecomp
      integer iuEdgelo0,iuEdgelo1
      integer iuEdgehi0,iuEdgehi1
      REAL_T uEdge(
     &           iuEdgelo0:iuEdgehi0,
     &           iuEdgelo1:iuEdgehi1,
     &           0:nuEdgecomp-1)
      integer ndivcomp
      integer idivlo0,idivlo1
      integer idivhi0,idivhi1
      REAL_T div(
     &           idivlo0:idivhi0,
     &           idivlo1:idivhi1,
     &           0:ndivcomp-1)
      integer igridIntlo0,igridIntlo1
      integer igridInthi0,igridInthi1
      REAL_T dx
      integer idir
      integer i0,i1
      integer c2fLo0,c2fLo1
      integer c2fHi0,c2fHi1
      integer comp
      REAL_T one_on_dx
      
      c2fLo0= 0*CHF_ID(0, idir)

      c2fLo1= 0*CHF_ID(1, idir)

      
      c2fHi0= 1*CHF_ID(0, idir)

      c2fHi1= 1*CHF_ID(1, idir)

      one_on_dx = one/dx
      do comp = 0, ndivcomp-1
         
      do i1 = igridIntlo1,igridInthi1
      do i0 = igridIntlo0,igridInthi0

            div(i0,i1, comp) = div(i0,i1, comp)
     &        + one_on_dx * 
     &        ( uEdge(i0 +c2fHi0,i1 +c2fHi1, comp)
     &        - uEdge(i0 +c2fLo0,i1 +c2fLo1, comp) )
         
      enddo
      enddo
      enddo
      return
      end
