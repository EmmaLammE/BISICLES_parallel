      subroutine CCDERIV(
     &           deriv
     &           ,iderivlo0,iderivlo1
     &           ,iderivhi0,iderivhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,iderivBoxlo0,iderivBoxlo1
     &           ,iderivBoxhi0,iderivBoxhi1
     &           ,dx
     &           ,derivDir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iderivlo0,iderivlo1
      integer iderivhi0,iderivhi1
      REAL*8 deriv(
     &           iderivlo0:iderivhi0,
     &           iderivlo1:iderivhi1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1)
      integer iderivBoxlo0,iderivBoxlo1
      integer iderivBoxhi0,iderivBoxhi1
      REAL*8 dx
      integer derivDir
      integer i0,i1
      integer ii0,ii1
      REAL*8 halfOnDx
      halfOnDx = (0.500d0)/dx
      ii0 = CHF_ID(0,derivDir)
      ii1 = CHF_ID(1,derivDir)
      do i1 = iderivBoxlo1,iderivBoxhi1
      do i0 = iderivBoxlo0,iderivBoxhi0
        deriv(i0,i1) = halfOnDx*(data(i0+ii0,i1+ii1) 
     &                                 -data(i0-ii0,i1-ii1) )
      enddo
      enddo
      return
      end
      subroutine CCDERIVMASK(
     &           deriv
     &           ,iderivlo0,iderivlo1
     &           ,iderivhi0,iderivhi1
     &           ,data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,thk
     &           ,ithklo0,ithklo1
     &           ,ithkhi0,ithkhi1
     &           ,iderivBoxlo0,iderivBoxlo1
     &           ,iderivBoxhi0,iderivBoxhi1
     &           ,dx
     &           ,derivDir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iderivlo0,iderivlo1
      integer iderivhi0,iderivhi1
      REAL*8 deriv(
     &           iderivlo0:iderivhi0,
     &           iderivlo1:iderivhi1)
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1)
      integer ithklo0,ithklo1
      integer ithkhi0,ithkhi1
      REAL*8 thk(
     &           ithklo0:ithkhi0,
     &           ithklo1:ithkhi1)
      integer iderivBoxlo0,iderivBoxlo1
      integer iderivBoxhi0,iderivBoxhi1
      REAL*8 dx
      integer derivDir
      integer i0,i1
      integer ii0,ii1
      REAL*8 halfOnDx, oneOnDx, thkl,thkr
      oneOnDx = (1.0d0)/dx
      halfOnDx = (0.500d0)/dx
      ii0 = CHF_ID(0,derivDir)
      ii1 = CHF_ID(1,derivDir)
      do i1 = iderivBoxlo1,iderivBoxhi1
      do i0 = iderivBoxlo0,iderivBoxhi0
      deriv(i0,i1) = (0.0d0)
      if (thk(i0,i1).gt.(0.0d0)) then
           thkl = thk(i0-ii0,i1-ii1)
           thkr = thk(i0+ii0,i1+ii1)
           if ((thkl.gt.(0.0d0)).and.(thkr.gt.(0.0d0))) then
              deriv(i0,i1) = halfOnDx*(data(i0+ii0,i1+ii1) 
     &             -data(i0-ii0,i1-ii1))
           else if (thkl.gt.(0.0d0)) then
              deriv(i0,i1) = oneOnDx * (data(i0,i1) 
     &             -data(i0-ii0,i1-ii1))
           else if (thkr.gt.(0.0d0)) then
              deriv(i0,i1) = oneOnDx *(data(i0+ii0,i1+ii1) 
     &             -data(i0,i1)) 
           end if
      end if
      enddo
      enddo
      return
      end
      subroutine FACEDERIV(
     &           deriv
     &           ,iderivlo0,iderivlo1
     &           ,iderivhi0,iderivhi1
     &           ,ccData
     &           ,iccDatalo0,iccDatalo1
     &           ,iccDatahi0,iccDatahi1
     &           ,iderivBoxlo0,iderivBoxlo1
     &           ,iderivBoxhi0,iderivBoxhi1
     &           ,dx
     &           ,derivDir
     &           ,faceDir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iderivlo0,iderivlo1
      integer iderivhi0,iderivhi1
      REAL*8 deriv(
     &           iderivlo0:iderivhi0,
     &           iderivlo1:iderivhi1)
      integer iccDatalo0,iccDatalo1
      integer iccDatahi0,iccDatahi1
      REAL*8 ccData(
     &           iccDatalo0:iccDatahi0,
     &           iccDatalo1:iccDatahi1)
      integer iderivBoxlo0,iderivBoxlo1
      integer iderivBoxhi0,iderivBoxhi1
      REAL*8 dx
      integer derivDir
      integer faceDir
      REAL*8 oneOnDx
      integer i0,i1
      integer ii0,ii1
      integer jj0,jj1
      oneOnDx = (1.0d0)/dx
      if (derivDir.eq.faceDir) then
         ii0 = CHF_ID(0,derivDir)
         ii1 = CHF_ID(1,derivDir)
      do i1 = iderivBoxlo1,iderivBoxhi1
      do i0 = iderivBoxlo0,iderivBoxhi0
            deriv(i0,i1) = oneOnDx*(ccData(i0,i1) 
     &                                     -ccData(i0-ii0,i1-ii1) )
      enddo
      enddo
      else
         ii0 = CHF_ID(0,derivDir)
         ii1 = CHF_ID(1,derivDir)
         jj0 = CHF_ID(0,faceDir)
         jj1 = CHF_ID(1,faceDir)
      do i1 = iderivBoxlo1,iderivBoxhi1
      do i0 = iderivBoxlo0,iderivBoxhi0
            deriv(i0,i1) = (0.250d0)*oneOnDx*(ccData(i0+ii0,i1+ii1)
     &                                    -ccData(i0-ii0,i1-ii1) 
     &                    +ccData(i0+ii0-jj0,i1+ii1-jj1)
     &                    -ccData(i0-ii0-jj0,i1-ii1-jj1))
      enddo
      enddo
       endif
      return
      end
