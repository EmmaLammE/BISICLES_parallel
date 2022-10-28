      subroutine MASKDOTPROD(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,y
     &           ,iylo0,iylo1
     &           ,iyhi0,iyhi1
     &           ,nycomp
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,nmaskcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,dotProd
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     &           iylo0:iyhi0,
     &           iylo1:iyhi1,
     &           0:nycomp-1)
      integer nmaskcomp
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      REAL*8 mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           0:nmaskcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 dotProd
      integer i0,i1
      integer nv, ncomp
      ncomp = nxcomp
      if (ncomp .ne. nycomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      dotProd = dotProd + x(i0,i1,nv) * y(i0,i1,nv) *
     &          mask(i0,i1,0)
      enddo
      enddo
      end do
      return
      end
      subroutine ARRAYPROD(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,y
     &           ,iylo0,iylo1
     &           ,iyhi0,iyhi1
     &           ,nycomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     &           iylo0:iyhi0,
     &           iylo1:iyhi1,
     &           0:nycomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nycomp .or. ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      z(i0,i1,nv) = x(i0,i1,nv) * y(i0,i1,nv)
      enddo
      enddo
      end do
      return
      end
      subroutine ARRAYDIV(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,y
     &           ,iylo0,iylo1
     &           ,iyhi0,iyhi1
     &           ,nycomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     &           iylo0:iyhi0,
     &           iylo1:iyhi1,
     &           0:nycomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nycomp .or. ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      z(i0,i1,nv) = x(i0,i1,nv) / y(i0,i1,nv)
      enddo
      enddo
      end do
      return
      end
      subroutine ARRAYSCL(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,c
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 c
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      z(i0,i1,nv) = c * x(i0,i1,nv)
      enddo
      enddo
      end do
      return
      end
       subroutine ARRAYABS(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      z(i0,i1,nv) = abs(x(i0,i1,nv))
      enddo
      enddo
      end do
      return
      end
      subroutine ADDCONST(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,b
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 b
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      z(i0,i1,nv) = x(i0,i1,nv) + b
      enddo
      enddo
      end do
      return
      end
      subroutine MASKSUMXP(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,nmaskcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,sumxp
     &           ,p
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nmaskcomp
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      REAL*8 mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           0:nmaskcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 sumxp
      integer p
      integer i0,i1
      integer nv, ncomp
      REAL*8 xpmask
      ncomp = nxcomp
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      xpmask = mask(i0,i1,0)*abs(x(i0,i1,nv))**p
      sumxp = sumxp+xpmask
      enddo
      enddo
      end do
      return
      end
      subroutine MASKMAXNORM(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,nmaskcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,m
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nmaskcomp
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      REAL*8 mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           0:nmaskcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 m
      integer i0,i1
      integer nv, ncomp
      REAL*8 absxmask
      ncomp = nxcomp
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      absxmask = abs(x(i0,i1,nv))*mask(i0,i1,0)
      m = max(abs(m), absxmask)
      enddo
      enddo
      end do
      return
      end
      subroutine MASKMIN(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,nmaskcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,m
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nmaskcomp
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      REAL*8 mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           0:nmaskcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 m
      integer i0,i1
      integer nv, ncomp
      ncomp = nxcomp
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      if (mask(i0,i1,0) .ne. 1) then
        m = min(m, x(i0,i1,nv))
      endif
      enddo
      enddo
      end do
      return
      end
      subroutine MASKWTDSQ(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,w
     &           ,iwlo0,iwlo1
     &           ,iwhi0,iwhi1
     &           ,nwcomp
     &           ,mask
     &           ,imasklo0,imasklo1
     &           ,imaskhi0,imaskhi1
     &           ,nmaskcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,norm
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nwcomp
      integer iwlo0,iwlo1
      integer iwhi0,iwhi1
      REAL*8 w(
     &           iwlo0:iwhi0,
     &           iwlo1:iwhi1,
     &           0:nwcomp-1)
      integer nmaskcomp
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      REAL*8 mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           0:nmaskcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 norm
      integer i0,i1
      integer nv, ncomp
      ncomp = nxcomp
      if(ncomp .ne. nwcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      norm = norm + 
     &       (x(i0,i1,nv)*w(i0,i1,nv)) * 
     &       (x(i0,i1,nv)*w(i0,i1,nv)) *
     &       mask(i0,i1,0)
      enddo
      enddo
      end do
      return
      end
      subroutine WTDSIGNSQ(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,w
     &           ,iwlo0,iwlo1
     &           ,iwhi0,iwhi1
     &           ,nwcomp
     &           ,id
     &           ,iidlo0,iidlo1
     &           ,iidhi0,iidhi1
     &           ,nidcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,norm
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nwcomp
      integer iwlo0,iwlo1
      integer iwhi0,iwhi1
      REAL*8 w(
     &           iwlo0:iwhi0,
     &           iwlo1:iwhi1,
     &           0:nwcomp-1)
      integer nidcomp
      integer iidlo0,iidlo1
      integer iidhi0,iidhi1
      REAL*8 id(
     &           iidlo0:iidhi0,
     &           iidlo1:iidhi1,
     &           0:nidcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 norm
      integer i0,i1
      integer nv, ncomp
      ncomp = nxcomp
      if(ncomp .ne. nwcomp .or. ncomp .ne. nidcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      norm = norm + 
     & (x(i0,i1,nv)*w(i0,i1,nv)*id(i0,i1,nv)) * 
     & (x(i0,i1,nv)*w(i0,i1,nv)*id(i0,i1,nv))
      enddo
      enddo
      end do
      return
      end
      subroutine ARRAYCOMP(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,c
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 c
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      if (x(i0,i1,nv) .ge. c) then
      z(i0,i1,nv) = 1.0
      else
      z(i0,i1,nv) = 0.0;
      end if
      enddo
      enddo
      end do
      return
      end
      subroutine INVWCHK(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,z
     &           ,izlo0,izlo1
     &           ,izhi0,izhi1
     &           ,nzcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,nonzero
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nzcomp
      integer izlo0,izlo1
      integer izhi0,izhi1
      REAL*8 z(
     &           izlo0:izhi0,
     &           izlo1:izhi1,
     &           0:nzcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      integer nonzero
      integer i0,i1
      integer nv, ncomp
      ncomp = nzcomp
      if(ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      if (x(i0,i1,nv) .ne. 0.0) then
      z(i0,i1,nv) = 1./x(i0,i1,nv)
      else
      z(i0,i1,nv) = 1.7e38
      nonzero = 0
      end if
      enddo
      enddo
      end do
      return
      end
      subroutine CONSTRCHK(
     &           c
     &           ,iclo0,iclo1
     &           ,ichi0,ichi1
     &           ,nccomp
     &           ,x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,m
     &           ,imlo0,imlo1
     &           ,imhi0,imhi1
     &           ,nmcomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,allpassed
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nccomp
      integer iclo0,iclo1
      integer ichi0,ichi1
      REAL*8 c(
     &           iclo0:ichi0,
     &           iclo1:ichi1,
     &           0:nccomp-1)
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nmcomp
      integer imlo0,imlo1
      integer imhi0,imhi1
      REAL*8 m(
     &           imlo0:imhi0,
     &           imlo1:imhi1,
     &           0:nmcomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      integer allpassed
      integer i0,i1
      integer nv, ncomp
      ncomp = nmcomp
      if(ncomp .ne. nccomp .or. ncomp .ne. nxcomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      if (c(i0,i1,nv) == 2) then
      if (x(i0,i1,nv) > 0) then
      m(i0,i1,nv) = 0.0
      else
      m(i0,i1,nv) = 1.0
      allpassed = 0
      end if
      else if (c(i0,i1,nv) == 1) then
      if (x(i0,i1,nv) .ge. 0) then
      m(i0,i1,nv) = 0.0
      else
      m(i0,i1,nv) = 1.0
      allpassed = 0
      end if
      else if (c(i0,i1,nv) == 0) then
      if (x(i0,i1,nv) < 0) then
      m(i0,i1,nv) = 0.0
      else
      m(i0,i1,nv) = 1.0
      allpassed = 0
      end if
      end if
      enddo
      enddo
      end do
      return
      end
      subroutine MINQUOT(
     &           x
     &           ,ixlo0,ixlo1
     &           ,ixhi0,ixhi1
     &           ,nxcomp
     &           ,y
     &           ,iylo0,iylo1
     &           ,iyhi0,iyhi1
     &           ,nycomp
     &           ,ireglo0,ireglo1
     &           ,ireghi0,ireghi1
     &           ,q
     &           ,nmin
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nxcomp
      integer ixlo0,ixlo1
      integer ixhi0,ixhi1
      REAL*8 x(
     &           ixlo0:ixhi0,
     &           ixlo1:ixhi1,
     &           0:nxcomp-1)
      integer nycomp
      integer iylo0,iylo1
      integer iyhi0,iyhi1
      REAL*8 y(
     &           iylo0:iyhi0,
     &           iylo1:iyhi1,
     &           0:nycomp-1)
      integer ireglo0,ireglo1
      integer ireghi0,ireghi1
      REAL*8 q
      integer nmin
      integer i0,i1
      integer nv, ncomp
      REAL*8 frac
      ncomp = nxcomp
      if(ncomp .ne. nycomp) then
      call MAYDAY_ERROR(
     &           )
      endif
      nmin = 1
      do nv = 0, ncomp-1
      do i1 = ireglo1,ireghi1
      do i0 = ireglo0,ireghi0
      if (y(i0,i1,nv) .ne. 0) then
      frac = x(i0,i1,nv) / y(i0,i1,nv)
      if (frac == q) then
      nmin = nmin + 1
      else if (frac < q) then
      q = frac
      nmin = 1
      end if
      end if
      enddo
      enddo
      end do
      return
      end
