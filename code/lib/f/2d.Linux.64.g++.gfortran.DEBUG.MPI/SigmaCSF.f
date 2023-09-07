      subroutine CENTEREDDIFF(
     & deltaPhi
     & ,ideltaPhilo0,ideltaPhilo1
     & ,ideltaPhihi0,ideltaPhihi1
     & ,phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,iderivBoxlo0,iderivBoxlo1
     & ,iderivBoxhi0,iderivBoxhi1
     & ,dx
     & ,dir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ideltaPhilo0,ideltaPhilo1
      integer ideltaPhihi0,ideltaPhihi1
      REAL*8 deltaPhi(
     & ideltaPhilo0:ideltaPhihi0,
     & ideltaPhilo1:ideltaPhihi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1)
      integer iderivBoxlo0,iderivBoxlo1
      integer iderivBoxhi0,iderivBoxhi1
      REAL*8 dx
      integer dir
      integer i0,i1
      integer ii0,ii1
      ii0 = CHF_ID(dir,0)
      ii1 = CHF_ID(dir,1)
      do i1 = iderivBoxlo1,iderivBoxhi1
      do i0 = iderivBoxlo0,iderivBoxhi0
        deltaPhi(i0,i1) = (phi(i0+ii0,i1+ii1)
     & -phi(i0-ii0,i1-ii1) )
      enddo
      enddo
      return
      end
      subroutine DEFINECELLGEOM(
     & deltaFactors
     & ,ideltaFactorslo0,ideltaFactorslo1
     & ,ideltaFactorshi0,ideltaFactorshi1
     & ,deltaH
     & ,ideltaHlo0,ideltaHlo1
     & ,ideltaHhi0,ideltaHhi1
     & ,deltaZb
     & ,ideltaZblo0,ideltaZblo1
     & ,ideltaZbhi0,ideltaZbhi1
     & ,dx
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,derivDir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ideltaFactorslo0,ideltaFactorslo1
      integer ideltaFactorshi0,ideltaFactorshi1
      REAL*8 deltaFactors(
     & ideltaFactorslo0:ideltaFactorshi0,
     & ideltaFactorslo1:ideltaFactorshi1)
      integer ideltaHlo0,ideltaHlo1
      integer ideltaHhi0,ideltaHhi1
      REAL*8 deltaH(
     & ideltaHlo0:ideltaHhi0,
     & ideltaHlo1:ideltaHhi1)
      integer ideltaZblo0,ideltaZblo1
      integer ideltaZbhi0,ideltaZbhi1
      REAL*8 deltaZb(
     & ideltaZblo0:ideltaZbhi0,
     & ideltaZblo1:ideltaZbhi1)
      REAL*8 dx(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer derivDir
      integer Hindex
      integer Zindex
      integer i0,i1
      REAL*8 sigma
      Hindex = ideltaHlo0
      Zindex = ideltaZblo0
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
        sigma = dx(0)*(i0 + (0.500d0))
        deltaFactors(i0,i1) = deltaH(Hindex,i1)
     & + deltaZb(Zindex,i1)
     & - sigma*deltaH(Hindex,i1)
      enddo
      enddo
      return
      end
      subroutine FIXFACEH(
     & xFaceH
     & ,ixFaceHlo0,ixFaceHlo1
     & ,ixFaceHhi0,ixFaceHhi1
     & ,yFaceH
     & ,iyFaceHlo0,iyFaceHlo1
     & ,iyFaceHhi0,iyFaceHhi1
     & ,cellH
     & ,icellHlo0,icellHlo1
     & ,icellHhi0,icellHhi1
     & ,igridboxlo0,igridboxlo1
     & ,igridboxhi0,igridboxhi1
     & ,xDir
     & ,yDir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ixFaceHlo0,ixFaceHlo1
      integer ixFaceHhi0,ixFaceHhi1
      REAL*8 xFaceH(
     & ixFaceHlo0:ixFaceHhi0,
     & ixFaceHlo1:ixFaceHhi1)
      integer iyFaceHlo0,iyFaceHlo1
      integer iyFaceHhi0,iyFaceHhi1
      REAL*8 yFaceH(
     & iyFaceHlo0:iyFaceHhi0,
     & iyFaceHlo1:iyFaceHhi1)
      integer icellHlo0,icellHlo1
      integer icellHhi0,icellHhi1
      REAL*8 cellH(
     & icellHlo0:icellHhi0,
     & icellHlo1:icellHhi1)
      integer igridboxlo0,igridboxlo1
      integer igridboxhi0,igridboxhi1
      integer xDir
      integer yDir
      integer i0,i1
      integer xFaceOff0,xFaceOff1
      integer yFaceOff0,yFaceOff1
      REAL*8 zeroVal;
      zeroVal = 1.0e-18
      xFaceOff0 = CHF_ID(xDir,0)
      xFaceOff1 = CHF_ID(xDir,1)
      yFaceOff0 = CHF_ID(yDir,0)
      yFaceOff1 = CHF_ID(yDir,1)
      do i1 = igridboxlo1,igridboxhi1
      do i0 = igridboxlo0,igridboxhi0
         if (abs(cellH(i0,i1)) .lt. zeroVal) then
            xFaceH(i0,i1) = (0.0d0)
            xFaceH(i0+xFaceOff0,i1+xFaceOff1) = (0.0d0)
            yFaceH(i0,i1) = (0.0d0)
            yFaceH(i0+yFaceOff0,i1+yFaceOff1) = (0.0d0)
         endif
      enddo
      enddo
      return
      end
      subroutine SURFACEHEIGHT(
     & zSurface
     & ,izSurfacelo0,izSurfacelo1
     & ,izSurfacehi0,izSurfacehi1
     & ,H
     & ,iHlo0,iHlo1
     & ,iHhi0,iHhi1
     & ,baseHeight
     & ,ibaseHeightlo0,ibaseHeightlo1
     & ,ibaseHeighthi0,ibaseHeighthi1
     & ,iceDensity
     & ,waterDensity
     & ,seaLevel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer izSurfacelo0,izSurfacelo1
      integer izSurfacehi0,izSurfacehi1
      REAL*8 zSurface(
     & izSurfacelo0:izSurfacehi0,
     & izSurfacelo1:izSurfacehi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     & iHlo0:iHhi0,
     & iHlo1:iHhi1)
      integer ibaseHeightlo0,ibaseHeightlo1
      integer ibaseHeighthi0,ibaseHeighthi1
      REAL*8 baseHeight(
     & ibaseHeightlo0:ibaseHeighthi0,
     & ibaseHeightlo1:ibaseHeighthi1)
      REAL*8 iceDensity
      REAL*8 waterDensity
      REAL*8 seaLevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 ratio, zs
      REAL*8 groundedHeight, floatingHeight
      ratio = (1.0d0) - (iceDensity/waterDensity)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         groundedHeight = H(i0,i1) + baseHeight(i0,i1)
         floatingHeight = seaLevel + ratio*H(i0,i1)
         zs = max(groundedHeight, floatingHeight);
         zSurface(i0,i1) = zs
      enddo
      enddo
      return
      end
      subroutine SETFLOATINGMASK(
     & floatingMask
     & ,ifloatingMasklo0,ifloatingMasklo1
     & ,ifloatingMaskhi0,ifloatingMaskhi1
     & ,Zsurf
     & ,iZsurflo0,iZsurflo1
     & ,iZsurfhi0,iZsurfhi1
     & ,zBase
     & ,izBaselo0,izBaselo1
     & ,izBasehi0,izBasehi1
     & ,H
     & ,iHlo0,iHlo1
     & ,iHhi0,iHhi1
     & ,anyFloating
     & ,iceDensity
     & ,waterDensity
     & ,seaLevel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifloatingMasklo0,ifloatingMasklo1
      integer ifloatingMaskhi0,ifloatingMaskhi1
      integer floatingMask(
     & ifloatingMasklo0:ifloatingMaskhi0,
     & ifloatingMasklo1:ifloatingMaskhi1)
      integer iZsurflo0,iZsurflo1
      integer iZsurfhi0,iZsurfhi1
      REAL*8 Zsurf(
     & iZsurflo0:iZsurfhi0,
     & iZsurflo1:iZsurfhi1)
      integer izBaselo0,izBaselo1
      integer izBasehi0,izBasehi1
      REAL*8 zBase(
     & izBaselo0:izBasehi0,
     & izBaselo1:izBasehi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     & iHlo0:iHhi0,
     & iHlo1:iHhi1)
      integer anyFloating
      REAL*8 iceDensity
      REAL*8 waterDensity
      REAL*8 seaLevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 eps
      eps = 1.0e-10
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (H(i0,i1).lt.eps) then
         if (zBase(i0,i1).lt.seaLevel) then
            anyFloating = 1;
            floatingMask(i0,i1) = (4)
         else
            floatingMask(i0,i1) = (8)
         end if
      else if (Zsurf(i0,i1) .gt. (zBase(i0,i1)
     & +H(i0,i1)+eps*H(i0,i1))) then
         anyFloating = 1;
         floatingMask(i0,i1) = (2)
      else
         floatingMask(i0,i1) = (1)
      endif
      enddo
      enddo
      return
      end
      subroutine PRESERVEMASK(
     & lsrf
     & ,ilsrflo0,ilsrflo1
     & ,ilsrfhi0,ilsrfhi1
     & ,usrf
     & ,iusrflo0,iusrflo1
     & ,iusrfhi0,iusrfhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,iceDensity
     & ,waterDensity
     & ,seaLevel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilsrflo0,ilsrflo1
      integer ilsrfhi0,ilsrfhi1
      REAL*8 lsrf(
     & ilsrflo0:ilsrfhi0,
     & ilsrflo1:ilsrfhi1)
      integer iusrflo0,iusrflo1
      integer iusrfhi0,iusrfhi1
      REAL*8 usrf(
     & iusrflo0:iusrfhi0,
     & iusrflo1:iusrfhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      REAL*8 iceDensity
      REAL*8 waterDensity
      REAL*8 seaLevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 ratio, oneOnRatio, lapa, lapb
      ratio = (iceDensity / waterDensity)
      oneOnRatio = (1.0d0)/ratio;
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (mask(i0,i1) .eq. (1)) then
         thck(i0,i1) =
     & max(usrf(i0,i1)-lsrf(i0,i1),
     & (1.0d-2))
         thck(i0,i1) = max(thck(i0,i1),
     & (seaLevel-topg(i0,i1))*oneOnRatio)
         lsrf(i0,i1) = topg(i0,i1);
         usrf(i0,i1) = lsrf(i0,i1) + thck(i0,i1)
      else if (mask(i0,i1) .eq. (8)) then
         topg(i0,i1) = max(topg(i0,i1),
     & seaLevel+(1.0d-2))
         lsrf(i0,i1) = topg(i0,i1)
         usrf(i0,i1) = topg(i0,i1)
         thck(i0,i1) = (0.0d0)
      else if (mask(i0,i1) .eq. (4)) then
         lsrf(i0,i1) = seaLevel;
         usrf(i0,i1) = seaLevel;
         topg(i0,i1) = min(topg(i0,i1),
     & seaLevel-(1.0d-2))
         thck(i0,i1) = (0.0d0)
      else
         topg(i0,i1) = min(
     & topg(i0,i1),
     & seaLevel - (2.0d0)*(1.0d-2))
         thck(i0,i1) =max(
     & usrf(i0,i1)-lsrf(i0,i1),
     & (1.0d-2))
         thck(i0,i1) = min(
     & thck(i0,i1),
     & (seaLevel-topg(i0,i1))*oneOnRatio-(1.0d-2))
         usrf(i0,i1) = (1-ratio)*thck(i0,i1) + seaLevel;
         lsrf(i0,i1) = usrf(i0,i1) - thck(i0,i1)
         topg(i0,i1) = min(topg(i0,i1),
     & lsrf(i0,i1) - (2.0d0)*(1.0d-2))
      end if
      enddo
      enddo
      return
      end
      subroutine PRESERVEMASKR(
     & lsrf
     & ,ilsrflo0,ilsrflo1
     & ,ilsrfhi0,ilsrfhi1
     & ,usrf
     & ,iusrflo0,iusrflo1
     & ,iusrfhi0,iusrfhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,iceDensity
     & ,waterDensity
     & ,seaLevel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilsrflo0,ilsrflo1
      integer ilsrfhi0,ilsrfhi1
      REAL*8 lsrf(
     & ilsrflo0:ilsrfhi0,
     & ilsrflo1:ilsrfhi1)
      integer iusrflo0,iusrflo1
      integer iusrfhi0,iusrfhi1
      REAL*8 usrf(
     & iusrflo0:iusrfhi0,
     & iusrflo1:iusrfhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      REAL*8 iceDensity
      REAL*8 waterDensity
      REAL*8 seaLevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 ratio, thck, oneOnRatio
      ratio = (iceDensity / waterDensity)
      oneOnRatio = (1.0d0)/ratio;
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (mask(i0,i1) .eq. (1)) then
         thck = max(usrf(i0,i1)-lsrf(i0,i1),
     & (1.0d-2))
         thck = max(thck, (seaLevel-topg(i0,i1))*oneOnRatio)
         lsrf(i0,i1) = topg(i0,i1);
         usrf(i0,i1) = lsrf(i0,i1) + thck;
      else if (mask(i0,i1) .eq. (8)) then
         topg(i0,i1) = max(topg(i0,i1),
     & seaLevel+(1.0d-2))
         lsrf(i0,i1) = topg(i0,i1);
         usrf(i0,i1) = topg(i0,i1);
      else if (mask(i0,i1) .eq. (4)) then
         lsrf(i0,i1) = seaLevel;
         usrf(i0,i1) = seaLevel;
         topg(i0,i1) = min(topg(i0,i1),
     & seaLevel-(1.0d-2))
      else
         topg(i0,i1) = min(topg(i0,i1),
     & seaLevel - (2.0d0)*(1.0d-2))
         thck = max(usrf(i0,i1)-lsrf(i0,i1),
     & (1.0d-2))
         thck = min(thck, (seaLevel-topg(i0,i1))*oneOnRatio
     & - (1.0d-2))
         usrf(i0,i1) = (1-ratio)*thck + seaLevel;
         lsrf(i0,i1) = usrf(i0,i1) - thck
         topg(i0,i1) = min(topg(i0,i1),
     & lsrf(i0,i1) - (2.0d0)*(1.0d-2))
      end if
      enddo
      enddo
      return
      end
      subroutine PRESERVEMASKO(
     & lsrf
     & ,ilsrflo0,ilsrflo1
     & ,ilsrfhi0,ilsrfhi1
     & ,usrf
     & ,iusrflo0,iusrflo1
     & ,iusrfhi0,iusrfhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,iceDensity
     & ,waterDensity
     & ,seaLevel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilsrflo0,ilsrflo1
      integer ilsrfhi0,ilsrfhi1
      REAL*8 lsrf(
     & ilsrflo0:ilsrfhi0,
     & ilsrflo1:ilsrfhi1)
      integer iusrflo0,iusrflo1
      integer iusrfhi0,iusrfhi1
      REAL*8 usrf(
     & iusrflo0:iusrfhi0,
     & iusrflo1:iusrfhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      REAL*8 iceDensity
      REAL*8 waterDensity
      REAL*8 seaLevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 oneOnRatio,usc,sslc
      REAL*8 tempThck
      oneOnRatio = (1.0d0) / ((1.0d0) - iceDensity / waterDensity)
      usc = (1.0d0) - oneOnRatio
      sslc = seaLevel * oneOnRatio
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         tempThck = max(usrf(i0,i1) - lsrf(i0,i1),(0.0d0))
         if (mask(i0,i1) .eq. (1)
     & .or. mask(i0,i1) .eq. (8)) then
         lsrf(i0,i1) = topg(i0,i1)
         if (usrf(i0,i1) .le. lsrf(i0,i1) ) then
            usrf(i0,i1) = lsrf(i0,i1) + tempThck
         endif
      else if (mask(i0,i1).eq.(4)) then
         if (usrf(i0,i1) .lt. topg(i0,i1)) then
            usrf(i0,i1) = topg(i0,i1)
         endif
         lsrf(i0,i1) = usrf(i0,i1)
      else
         lsrf(i0,i1) = max(topg(i0,i1) + (1.0d-2),
     & usrf(i0,i1) * usc - sslc)
         if (usrf(i0,i1) .le. lsrf(i0,i1) ) then
            usrf(i0,i1) = lsrf(i0,i1) + tempThck
         endif
      end if
      enddo
      enddo
      return
      end
      subroutine ODESTROYSCFI(
     & thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,count
     & ,xdir
     & ,ydir
     & ,thresh
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer count
      integer xdir
      integer ydir
      REAL*8 thresh
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ix0,ix1
      integer iy0,iy1
      ix0 = CHF_ID(xdir,0)
                ix1 = CHF_ID(xdir,1)
      iy0 = CHF_ID(ydir,0)
                iy1 = CHF_ID(ydir,1)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if ( mask(i0,i1).eq.(2) ) then
         if ( (thck(i0-ix0,i1-ix1).le.thresh)
     & .and.(thck(i0+ix0,i1+ix1).le.thresh)
     & .and.(thck(i0-iy0,i1-iy1).le.thresh)
     & .and.(thck(i0+iy0,i1+iy1).le.thresh)) then
            count = count + 1
            thck(i0,i1)=0.0;
            mask(i0,i1)=(4)
         end if
      end if
      enddo
      enddo
      end subroutine
      subroutine DESTROYSCFI(
     & thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,count
     & ,dir
     & ,thresh
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer count
      integer dir
      REAL*8 thresh
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer ii0,ii1
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if ( mask(i0,i1).eq.(2) ) then
         if ((thck(i0-ii0,i1-ii1).le.thresh)
     & .and.(thck(i0+ii0,i1+ii1).le.thresh)) then
            count = count + 1
            thck(i0,i1)=0.0;
            mask(i0,i1)=(4)
         end if
      end if
      enddo
      enddo
      end subroutine
      subroutine SSETOPENSURFACE(
     & surf
     & ,isurflo0,isurflo1
     & ,isurfhi0,isurfhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,seaLevel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isurflo0,isurflo1
      integer isurfhi0,isurfhi1
      REAL*8 surf(
     & isurflo0:isurfhi0,
     & isurflo1:isurfhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      REAL*8 seaLevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if ( mask(i0,i1).eq. (4) ) then
         surf(i0,i1) = seaLevel
      else if( mask(i0,i1).eq. (8) ) then
         surf(i0,i1) = topg(i0,i1)
      end if
      enddo
      enddo
      return
      end
      subroutine SGLGRADS(
     & grads
     & ,igradslo0,igradslo1
     & ,igradshi0,igradshi1
     & ,surff
     & ,isurfflo0,isurfflo1
     & ,isurffhi0,isurffhi1
     & ,thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,surf
     & ,isurflo0,isurflo1
     & ,isurfhi0,isurfhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,ratio
     & ,seaLevel
     & ,dx
     & ,dir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,ifaceBoxlo0,ifaceBoxlo1
     & ,ifaceBoxhi0,ifaceBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradslo0,igradslo1
      integer igradshi0,igradshi1
      REAL*8 grads(
     & igradslo0:igradshi0,
     & igradslo1:igradshi1)
      integer isurfflo0,isurfflo1
      integer isurffhi0,isurffhi1
      REAL*8 surff(
     & isurfflo0:isurffhi0,
     & isurfflo1:isurffhi1)
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      integer isurflo0,isurflo1
      integer isurfhi0,isurfhi1
      REAL*8 surf(
     & isurflo0:isurfhi0,
     & isurflo1:isurfhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      REAL*8 ratio
      REAL*8 seaLevel
      REAL*8 dx
      integer dir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ifaceBoxlo0,ifaceBoxlo1
      integer ifaceBoxhi0,ifaceBoxhi1
      integer i0,i1
      integer ii0,ii1
      integer maskP, maskL, maskR
      integer LISH
      REAL*8 oneOnDx
      REAL*8 topgf, thckf
      LISH = IOR((1), (2))
      oneOnDx = (1.0d0)/dx
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = ifaceboxlo1,ifaceboxhi1
      do i0 = ifaceboxlo0,ifaceboxhi0
      maskL = mask(i0-ii0,i1-ii1)
      maskR = mask(i0,i1)
      if (1.eq.0 .and. IOR(maskL,maskR) .eq. LISH) then
         topgf = (0.500d0)*(topg( i0,i1)
     & + topg( i0-ii0,i1-ii1))
         surff(i0,i1) = seaLevel
     & + ((1.0d0)-ratio)*(seaLevel - topgf)/ratio
      else
         surff(i0,i1) = (0.500d0) *
     & (surf(i0-ii0,i1-ii1) + surf(i0,i1))
      end if
      if (maskL.eq.(8)) then
         surff(i0,i1) =
     & min(surff(i0,i1), surf(i0,i1))
      end if
      if (maskR.eq.(8)) then
         surff(i0,i1) =
     & min(surff(i0,i1), surf(i0-ii0,i1-ii1))
      end if
      enddo
      enddo
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (iand(mask(i0,i1),LISH).gt.0) then
         grads(i0,i1) = oneOnDx *
     & (- surff(i0,i1) + surff(i0+ii0,i1+ii1))
      else
         grads(i0,i1) = (0.0d0)
      end if
      enddo
      enddo
      return
      end
      subroutine THICKNESSOVERFLOTATION(
     & p
     & ,iplo0,iplo1
     & ,iphi0,iphi1
     & ,thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,rhoi
     & ,rhoo
     & ,sealevel
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
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      integer itopglo0,itopglo1
      integer itopghi0,itopghi1
      REAL*8 topg(
     & itopglo0:itopghi0,
     & itopglo1:itopghi1)
      REAL*8 rhoi
      REAL*8 rhoo
      REAL*8 sealevel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 ratio
      REAL*8 eps
      eps = 1.0e-10
      ratio = rhoo / rhoi
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      p(i0,i1) = (0.0d0)
      if (thck(i0,i1).ge.eps) then
          p(i0,i1) = min( thck(i0,i1),
     & max((0.0d0),
     & (topg(i0,i1)-sealevel) * ratio
     & + thck(i0,i1)))
      endif
      enddo
      enddo
      return
      end
