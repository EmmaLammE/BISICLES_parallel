#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine INCREMENTLAPLACIAN(
     &           lapPhi
     &           ,ilapPhilo0,ilapPhilo1
     &           ,ilapPhihi0,ilapPhihi1
     &           ,nlapPhicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,igridBoxlo0,igridBoxlo1
     &           ,igridBoxhi0,igridBoxhi1
     &           ,dir
     &           ,factor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1
      integer ilapPhihi0,ilapPhihi1
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           0:nlapPhicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer dir
      REAL_T factor
      integer i0,i1
      integer ii0,ii1
      integer comp
      REAL_T thisLap
      
      ii0=CHF_ID(0,dir)

      ii1=CHF_ID(1,dir)

      do comp=0, nphicomp-1
         
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

            thisLap = phi(i0 +ii0,i1 +ii1, comp)
     &              + phi(i0 -ii0,i1 -ii1, comp)
     &           -two*phi(i0,i1, comp)
            lapPhi(i0,i1, comp) =
     &           lapPhi(i0,i1, comp) + factor*thisLap
         
      enddo
      enddo
      enddo
      return
      end
      subroutine INCREMENTLOSIDELAPLACIAN(
     &           lapPhi
     &           ,ilapPhilo0,ilapPhilo1
     &           ,ilapPhihi0,ilapPhihi1
     &           ,nlapPhicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,igridBoxlo0,igridBoxlo1
     &           ,igridBoxhi0,igridBoxhi1
     &           ,dir
     &           ,factor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1
      integer ilapPhihi0,ilapPhihi1
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           0:nlapPhicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer dir
      REAL_T factor
      integer i0,i1, comp
      integer ii0,ii1
      REAL_T thisLap
      
      ii0=CHF_ID(0,dir)

      ii1=CHF_ID(1,dir)

      do comp=0, nphicomp-1
         
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

            thisLap =
     &        two*phi(i0,i1,comp)
     &        - five*phi(i0 +ii0,i1 +ii1,comp)
     &        + four*phi(i0 +2*ii0,i1 +2*ii1,comp)
     &        - phi(i0 +3*ii0,i1 +3*ii1,comp)
            lapPhi(i0,i1,comp) =
     &           lapPhi(i0,i1,comp) + factor*thisLap
         
      enddo
      enddo
      enddo
      return
      end
      subroutine INCREMENTHISIDELAPLACIAN(
     &           lapPhi
     &           ,ilapPhilo0,ilapPhilo1
     &           ,ilapPhihi0,ilapPhihi1
     &           ,nlapPhicomp
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,igridBoxlo0,igridBoxlo1
     &           ,igridBoxhi0,igridBoxhi1
     &           ,dir
     &           ,factor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlapPhicomp
      integer ilapPhilo0,ilapPhilo1
      integer ilapPhihi0,ilapPhihi1
      REAL_T lapPhi(
     &           ilapPhilo0:ilapPhihi0,
     &           ilapPhilo1:ilapPhihi1,
     &           0:nlapPhicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer dir
      REAL_T factor
      integer i0,i1, comp
      integer ii0,ii1
      REAL_T thisLap
      
      ii0=CHF_ID(0,dir)

      ii1=CHF_ID(1,dir)

      do comp=0, nphicomp-1
         
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0

            thisLap =
     &        two*phi(i0,i1,comp)
     &        - five*phi(i0 -ii0,i1 -ii1,comp)
     &        + four*phi(i0 -2*ii0,i1 -2*ii1,comp)
     &        - phi(i0 -3*ii0,i1 -3*ii1,comp)
            lapPhi(i0,i1,comp) =
     &           lapPhi(i0,i1,comp) + factor*thisLap
         
      enddo
      enddo
      enddo
      return
      end
