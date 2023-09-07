      subroutine UNDIVIDEDGRAD(
     & dq
     & ,idqlo0,idqlo1
     & ,idqhi0,idqhi1
     & ,q
     & ,iqlo0,iqlo1
     & ,iqhi0,iqhi1
     & ,nqcomp
     & ,idInteriorlo0,idInteriorlo1
     & ,idInteriorhi0,idInteriorhi1
     & ,iloedgelo0,iloedgelo1
     & ,iloedgehi0,iloedgehi1
     & ,ihiedgelo0,ihiedgelo1
     & ,ihiedgehi0,ihiedgehi1
     & ,idir
     & ,haslo
     & ,hashi
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idqlo0,idqlo1
      integer idqhi0,idqhi1
      REAL*8 dq(
     & idqlo0:idqhi0,
     & idqlo1:idqhi1)
      integer nqcomp
      integer iqlo0,iqlo1
      integer iqhi0,iqhi1
      REAL*8 q(
     & iqlo0:iqhi0,
     & iqlo1:iqhi1,
     & 0:nqcomp-1)
      integer idInteriorlo0,idInteriorlo1
      integer idInteriorhi0,idInteriorhi1
      integer iloedgelo0,iloedgelo1
      integer iloedgehi0,iloedgehi1
      integer ihiedgelo0,ihiedgelo1
      integer ihiedgehi0,ihiedgehi1
      integer idir
      integer haslo
      integer hashi
      integer ldir, i,j
      integer ioff,joff
      integer n, ncomp
      ncomp = nqcomp
      ioff = CHF_ID(0,idir)
      joff = CHF_ID(1,idir)
      do j = idInteriorlo1,idInteriorhi1
      do i = idInteriorlo0,idInteriorhi0
      dq(i,j) = 0
      do n = 0, ncomp-1
       dq(i,j) = max( abs(dq(i,j)),
     & (0.500d0)*abs( q(i+ioff,j+joff,n)
     & - q(i-ioff,j-joff,n) ))
      enddo
      enddo
      enddo
      if (haslo .eq. 1) then
      do j = iloedgelo1,iloedgehi1
      do i = iloedgelo0,iloedgehi0
         dq(i,j)=0
         do n = 0, ncomp-1
           dq(i,j) = max( abs(dq(i,j)),
     & abs( q(i+ioff,j+joff, n) - q(i,j, n)))
         enddo
      enddo
      enddo
      endif
      if (hashi .eq. 1) then
      do j = ihiedgelo1,ihiedgehi1
      do i = ihiedgelo0,ihiedgehi0
         dq(i,j)=0
         do n = 0, ncomp-1
           dq(i,j) = max( abs(dq(i,j)),
     & abs(q(i,j, n) -
     & q(i-ioff,j-joff, n)))
         enddo
      enddo
      enddo
      endif
      return
      end
      subroutine SETONMASK(
     & fab
     & ,ifablo0,ifablo1
     & ,ifabhi0,ifabhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,maskval
     & ,fabval
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     & ifablo0:ifabhi0,
     & ifablo1:ifabhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer maskval
      REAL*8 fabval
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (mask(i0,i1).eq.maskval) then
         fab(i0,i1) = fabval
      end if
      enddo
      enddo
      return
      end
      subroutine SWEEPCONNECTED2D(
     & fab
     & ,ifablo0,ifablo1
     & ,ifabhi0,ifabhi1
     & ,conn
     & ,iconnlo0,iconnlo1
     & ,iconnhi0,iconnhi1
     & ,tol
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     & ifablo0:ifabhi0,
     & ifablo1:ifabhi1)
      integer iconnlo0,iconnlo1
      integer iconnhi0,iconnhi1
      REAL*8 conn(
     & iconnlo0:iconnhi0,
     & iconnlo1:iconnhi1)
      REAL*8 tol
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i,j
      do j = iboxlo1,iboxhi1
         do i = iboxlo0,iboxhi0,1
            if ((conn(i,j).gt.tol).and.(conn(i-1,j).gt.tol)) then
               fab(i,j) = max(fab(i,j),fab(i-1,j) )
            end if
         end do
         do i = iboxhi0,iboxlo0,-1
            if ((conn(i,j).gt.tol).and.(conn(i+1,j).gt.tol)) then
               fab(i,j) = max(fab(i,j),fab(i+1,j) )
            end if
         end do
      end do
      do i = iboxlo0,iboxhi0
         do j = iboxlo1,iboxhi1,1
            if ((conn(i,j).gt.tol).and.(conn(i,j-1).gt.tol)) then
               fab(i,j) = max(fab(i,j),fab(i,j-1) )
            end if
         end do
         do j = iboxhi1,iboxlo1,-1
            if ((conn(i,j).gt.tol).and.(conn(i,j+1).gt.tol)) then
               fab(i,j) = max(fab(i,j),fab(i,j+1) )
            end if
         end do
      end do
      return
      end
      subroutine SETFLOATINGBETA(
     & beta
     & ,ibetalo0,ibetalo1
     & ,ibetahi0,ibetahi1
     & ,floatingMask
     & ,ifloatingMasklo0,ifloatingMasklo1
     & ,ifloatingMaskhi0,ifloatingMaskhi1
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ibetalo0,ibetalo1
      integer ibetahi0,ibetahi1
      REAL*8 beta(
     & ibetalo0:ibetahi0,
     & ibetalo1:ibetahi1)
      integer ifloatingMasklo0,ifloatingMasklo1
      integer ifloatingMaskhi0,ifloatingMaskhi1
      integer floatingMask(
     & ifloatingMasklo0:ifloatingMaskhi0,
     & ifloatingMasklo1:ifloatingMaskhi1)
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer i0,i1
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0
      if (floatingMask(i0,i1).eq.(2)) then
         beta(i0,i1) = (0.0d0)
      else if (floatingMask(i0,i1).eq.(4)) then
         beta(i0,i1) = max(1.0d+2,beta(i0,i1))
      else if (floatingMask(i0,i1).eq.(8)) then
         beta(i0,i1) = max(1.0d+2,beta(i0,i1))
      else
         beta(i0,i1) = max(1.0d-10,beta(i0,i1))
      end if
      enddo
      enddo
      return
      end
      subroutine CORNFORDCORRECTION(
     & grad
     & ,igradlo0,igradlo1
     & ,igradhi0,igradhi1
     & ,Zsurf
     & ,iZsurflo0,iZsurflo1
     & ,iZsurfhi0,iZsurfhi1
     & ,floatingMask
     & ,ifloatingMasklo0,ifloatingMasklo1
     & ,ifloatingMaskhi0,ifloatingMaskhi1
     & ,dir
     & ,dx
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradlo0,igradlo1
      integer igradhi0,igradhi1
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1)
      integer iZsurflo0,iZsurflo1
      integer iZsurfhi0,iZsurfhi1
      REAL*8 Zsurf(
     & iZsurflo0:iZsurfhi0,
     & iZsurflo1:iZsurfhi1)
      integer ifloatingMasklo0,ifloatingMasklo1
      integer ifloatingMaskhi0,ifloatingMaskhi1
      integer floatingMask(
     & ifloatingMasklo0:ifloatingMaskhi0,
     & ifloatingMasklo1:ifloatingMaskhi1)
      integer dir
      REAL*8 dx
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer i0,i1
      integer ii0,ii1
      integer nstep
      integer radius
      REAL*8 oneOnDx
      oneOnDx = (1.0d0)/dx
      radius = 1
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0
        if (floatingMask(i0,i1).eq.(2)) then
           if (floatingMask(i0+ii0,i1+ii1).eq.(1)) then
              do nstep=0, radius-1
                 if (((i0-nstep*ii0.ge.igradlo0).AND.(i1-nstep*ii1.ge.ig
     &radlo1))) then
                    grad(i0-nstep*ii0,i1-nstep*ii1) = oneOnDx*(zSurf(i0-
     &nstep*ii0,i1-nstep*ii1) - zSurf(i0-(nstep+1)*ii0,i1-(nstep+1)*ii1)
     &)
              endif
              enddo
              do nstep=1, radius
                 if (((i0+nstep*ii0.le.igradhi0).AND.(i1+nstep*ii1.le.ig
     &radhi1))) then
                    grad(i0+nstep*ii0,i1+nstep*ii1) = oneOnDx*(zSurf(i0+
     &(nstep+1)*ii0,i1+(nstep+1)*ii1) - zSurf(i0+(nstep)*ii0,i1+(nstep)*
     &ii1))
                 endif
              enddo
           else if (floatingMask(i0-ii0,i1-ii1).eq.(1)) then
              do nstep=1, radius
                 if (((i0-nstep*ii0.ge.igradlo0).AND.(i1-nstep*ii1.ge.ig
     &radlo1))) then
                    grad(i0-nstep*ii0,i1-nstep*ii1) = oneOnDx*(zSurf(i0-
     &nstep*ii0,i1-nstep*ii1) - zSurf(i0-(nstep+1)*ii0,i1-(nstep+1)*ii1)
     &)
                 endif
              enddo
              do nstep=0, radius-1
                 if (((i0+nstep*ii0.le.igradhi0).AND.(i1+nstep*ii1.le.ig
     &radhi1))) then
                    grad(i0+nstep*ii0,i1+nstep*ii1) = oneOnDx*(zSurf(i0+
     &(nstep+1)*ii0,i1+(nstep+1)*ii1) - zSurf(i0+(nstep)*ii0,i1+(nstep)*
     &ii1))
                 endif
              enddo
           endif
        endif
      enddo
      enddo
      return
      end
      subroutine SETOPENSURFACE(
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
      subroutine SETICEFREEVAL(
     & fab
     & ,ifablo0,ifablo1
     & ,ifabhi0,ifabhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,val
     & ,icellBoxlo0,icellBoxlo1
     & ,icellBoxhi0,icellBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifablo0,ifablo1
      integer ifabhi0,ifabhi1
      REAL*8 fab(
     & ifablo0:ifabhi0,
     & ifablo1:ifabhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      REAL*8 val
      integer icellBoxlo0,icellBoxlo1
      integer icellBoxhi0,icellBoxhi1
      integer i0,i1
      integer maskc
      do i1 = icellBoxlo1,icellBoxhi1
      do i0 = icellBoxlo0,icellBoxhi0
      maskc = mask(i0,i1)
      if ((maskc.eq.(4)).or.(maskc.eq.(8))) then
         fab(i0,i1) = val
      end if
      enddo
      enddo
      return
      end
      subroutine SETFRONTFACEVT(
     & facevt
     & ,ifacevtlo0,ifacevtlo1
     & ,ifacevthi0,ifacevthi1
     & ,thck
     & ,ithcklo0,ithcklo1
     & ,ithckhi0,ithckhi1
     & ,usrf
     & ,iusrflo0,iusrflo1
     & ,iusrfhi0,iusrfhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,k
     & ,factor
     & ,icellBoxlo0,icellBoxlo1
     & ,icellBoxhi0,icellBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacevtlo0,ifacevtlo1
      integer ifacevthi0,ifacevthi1
      REAL*8 facevt(
     & ifacevtlo0:ifacevthi0,
     & ifacevtlo1:ifacevthi1)
      integer ithcklo0,ithcklo1
      integer ithckhi0,ithckhi1
      REAL*8 thck(
     & ithcklo0:ithckhi0,
     & ithcklo1:ithckhi1)
      integer iusrflo0,iusrflo1
      integer iusrfhi0,iusrfhi1
      REAL*8 usrf(
     & iusrflo0:iusrfhi0,
     & iusrflo1:iusrfhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer k
      REAL*8 factor
      integer icellBoxlo0,icellBoxlo1
      integer icellBoxhi0,icellBoxhi1
      integer i0,i1
      integer ii0,ii1
      integer maskc,maskr,maskl
      logical icel,icer,icec
      ii0 = CHF_ID(k,0)
                ii1 = CHF_ID(k,1)
      do i1 = icellBoxlo1,icellBoxhi1
      do i0 = icellBoxlo0,icellBoxhi0
      maskc = mask(i0,i1)
      icec = ((maskc.eq.(1)).or.(maskc.eq.(2)))
      if (icec) then
         maskl = mask(i0-ii0,i1-ii1)
         maskr = mask(i0+ii0,i1+ii1)
         icel = ((maskl.eq.(1)).or.(maskl.eq.(2)))
         icer = ((maskr.eq.(1)).or.(maskr.eq.(2)))
         if ( icer .and. (.not.icel)) then
            facevt(i0,i1) = facevt(i0+ii0,i1+ii1)
     & - factor * thck(i0,i1)
     & * (usrf(i0+ii0,i1+ii1) - usrf(i0,i1))
         end if
         if ( icel .and. (.not.icer)) then
            facevt(i0+ii0,i1+ii1) = facevt(i0,i1)
     & + factor * thck(i0,i1)
     & * (usrf(i0,i1)-usrf(i0-ii0,i1-ii1))
         end if
      end if
      enddo
      enddo
      return
      end
      subroutine EXTRAPTOMARGIN(
     & uf
     & ,iuflo0,iuflo1
     & ,iufhi0,iufhi1
     & ,vface
     & ,ivfacelo0,ivfacelo1
     & ,ivfacehi0,ivfacehi1
     & ,ufin
     & ,iufinlo0,iufinlo1
     & ,iufinhi0,iufinhi1
     & ,uc
     & ,iuclo0,iuclo1
     & ,iuchi0,iuchi1
     & ,usrf
     & ,iusrflo0,iusrflo1
     & ,iusrfhi0,iusrfhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,thk
     & ,ithklo0,ithklo1
     & ,ithkhi0,ithkhi1
     & ,dir
     & ,ifaceboxlo0,ifaceboxlo1
     & ,ifaceboxhi0,ifaceboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iuflo0,iuflo1
      integer iufhi0,iufhi1
      REAL*8 uf(
     & iuflo0:iufhi0,
     & iuflo1:iufhi1)
      integer ivfacelo0,ivfacelo1
      integer ivfacehi0,ivfacehi1
      REAL*8 vface(
     & ivfacelo0:ivfacehi0,
     & ivfacelo1:ivfacehi1)
      integer iufinlo0,iufinlo1
      integer iufinhi0,iufinhi1
      REAL*8 ufin(
     & iufinlo0:iufinhi0,
     & iufinlo1:iufinhi1)
      integer iuclo0,iuclo1
      integer iuchi0,iuchi1
      REAL*8 uc(
     & iuclo0:iuchi0,
     & iuclo1:iuchi1)
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
      integer ithklo0,ithklo1
      integer ithkhi0,ithkhi1
      REAL*8 thk(
     & ithklo0:ithkhi0,
     & ithklo1:ithkhi1)
      integer dir
      integer ifaceboxlo0,ifaceboxlo1
      integer ifaceboxhi0,ifaceboxhi1
      REAL*8 f,hl,hr
      integer i0,i1
      integer ii0,ii1
      logical icel,icer
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = ifaceboxlo1,ifaceboxhi1
      do i0 = ifaceboxlo0,ifaceboxhi0
      hl = thk(i0-ii0,i1-ii1)
      hr = thk(i0,i1)
      icel = (hl.gt.(1.0d-2))
      icer = (hr.gt.(1.0d-2))
      vface(i0,i1) = (1.0d0)
      uf(i0,i1) = ufin(i0,i1)
      if (icel.and.(.not.icer)) then
         uf(i0,i1) = max((2.0d0) * uc(i0-ii0,i1-ii1)
     & - ufin(i0-ii0,i1-ii1), ufin(i0,i1))
         vface(i0,i1) = (1.0d0)/hl * max((0.0d0),
     & usrf(i0-ii0,i1-ii1) - topg(i0,i1))
      else if (icer.and.(.not.icel)) then
         uf(i0,i1) = min((2.0d0) * uc(i0,i1)
     & - ufin(i0+ii0,i1+ii1),ufin(i0,i1))
         vface(i0,i1) = (1.0d0)/hr * max((0.0d0),
     & usrf(i0,i1) - topg(i0-ii0,i1-ii1))
      end if
      if (vface(i0,i1) .lt. (1.0d0) ) then
         uf(i0,i1) = uf(i0,i1) * vface(i0,i1)
      end if
      enddo
      enddo
      return
      end
      subroutine GLCORRECTION(
     & rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,H
     & ,iHlo0,iHlo1
     & ,iHhi0,iHhi1
     & ,Zsurf
     & ,iZsurflo0,iZsurflo1
     & ,iZsurfhi0,iZsurfhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,dir
     & ,dx
     & ,rhog
     & ,igridBoxlo0,igridBoxlo1
     & ,igridBoxhi0,igridBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     & iHlo0:iHhi0,
     & iHlo1:iHhi1)
      integer iZsurflo0,iZsurflo1
      integer iZsurfhi0,iZsurfhi1
      REAL*8 Zsurf(
     & iZsurflo0:iZsurfhi0,
     & iZsurflo1:iZsurfhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      integer dir
      REAL*8 dx
      REAL*8 rhog
      integer igridBoxlo0,igridBoxlo1
      integer igridBoxhi0,igridBoxhi1
      integer i0,i1
      integer ii0,ii1
      integer lmask,cmask,rmask
      integer LI,SH,OS,OL
      REAL*8 halfOnDx, sl,sr,sc,hl,hr,hc,tmprhs
      SH = (2)
      LI = (1)
      OS = (4)
      OL = (8)
      halfOnDx = (0.500d0)/dx * rhog
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = igridBoxlo1,igridBoxhi1
      do i0 = igridBoxlo0,igridBoxhi0
      cmask = mask(i0,i1)
      rmask = mask(i0+ii0,i1+ii1)
      lmask = mask(i0-ii0,i1-ii1)
      sc = zSurf(i0,i1)
      sl = zSurf(i0-ii0,i1-ii1)
      sr = zSurf(i0+ii0,i1+ii1)
      hc = H(i0,i1)
      hl = H(i0-ii0,i1-ii1)
      hr = H(i0+ii0,i1+ii1)
      if (cmask.eq.SH) then
         if (lmask.eq.LI) then
            if ( (rmask.eq.SH).or.(rmask.eq.OS) ) then
               tmprhs = halfOnDx * (sr-sc) * (hr+hc)
                  rhs(i0,i1) = tmprhs
            else
               rhs(i0,i1) = (0.0d0)
            end if
         else if (rmask.eq.LI) then
            if ( (lmask.eq.SH).or.(lmask.eq.OS) ) then
               tmprhs = halfOnDx * (sc-sl) * (hc+hl)
                  rhs(i0,i1) = tmprhs
            else
               rhs(i0,i1) = (0.0d0)
            endif
         end if
      else if (cmask.eq.LI) then
         if (lmask.eq.SH) then
            if (rmask.eq.LI) then
               tmprhs = halfOnDx * (sr-sc) * (hr+hc)
               if ( (tmprhs.gt.(0.0d0)) .and.
     & (tmprhs.gt.rhs(i0,i1))) then
                  rhs(i0,i1) = tmprhs
               end if
            else
               rhs(i0,i1) = rhs(i0,i1)
            end if
         else if (rmask.eq.SH) then
            if (lmask.eq.LI) then
               tmprhs = halfOnDx * (sc-sl) * (hc+hl)
               if ( (tmprhs.lt.(0.0d0)) .and.
     & (tmprhs.lt.rhs(i0,i1))) then
                  rhs(i0,i1) = tmprhs
               end if
            else
               rhs(i0,i1) = rhs(i0,i1)
            endif
         end if
      end if
      enddo
      enddo
      return
      end
      subroutine SETTHICKDIFF(
     & D
     & ,iDlo0,iDlo1
     & ,iDhi0,iDhi1
     & ,C
     & ,iClo0,iClo1
     & ,iChi0,iChi1
     & ,H
     & ,iHlo0,iHlo1
     & ,iHhi0,iHhi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,rg
     & ,dir
     & ,ifaceBoxlo0,ifaceBoxlo1
     & ,ifaceBoxhi0,ifaceBoxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iDlo0,iDlo1
      integer iDhi0,iDhi1
      REAL*8 D(
     & iDlo0:iDhi0,
     & iDlo1:iDhi1)
      integer iClo0,iClo1
      integer iChi0,iChi1
      REAL*8 C(
     & iClo0:iChi0,
     & iClo1:iChi1)
      integer iHlo0,iHlo1
      integer iHhi0,iHhi1
      REAL*8 H(
     & iHlo0:iHhi0,
     & iHlo1:iHhi1)
      integer imasklo0,imasklo1
      integer imaskhi0,imaskhi1
      integer mask(
     & imasklo0:imaskhi0,
     & imasklo1:imaskhi1)
      REAL*8 rg
      integer dir
      integer ifaceBoxlo0,ifaceBoxlo1
      integer ifaceBoxhi0,ifaceBoxhi1
      integer i0,i1
      integer ii0,ii1
      REAL*8 oneOnD
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = ifaceBoxlo1,ifaceBoxhi1
      do i0 = ifaceBoxlo0,ifaceBoxhi0
      if ( (mask(i0,i1).eq.(1)) .or.
     & (mask(i0-ii0,i1-ii1).eq.(1)) ) then
      oneOnD = 0.5 * ( C(i0,i1) / H(i0,i1)**2
     & + C (i0-ii0,i1-ii1) / H (i0-ii0,i1-ii1)**2)
         D(i0,i1) = rg / oneOnD
      else
         D(i0,i1) = 0.0;
      end if
      enddo
      enddo
      return
      end
      subroutine SUBTRACTDVEL(
     & faceVel
     & ,ifaceVello0,ifaceVello1
     & ,ifaceVelhi0,ifaceVelhi1
     & ,cellH
     & ,icellHlo0,icellHlo1
     & ,icellHhi0,icellHhi1
     & ,faceH
     & ,ifaceHlo0,ifaceHlo1
     & ,ifaceHhi0,ifaceHhi1
     & ,faceD
     & ,ifaceDlo0,ifaceDlo1
     & ,ifaceDhi0,ifaceDhi1
     & ,dir
     & ,dx
     & ,ifaceBoxlo0,ifaceBoxlo1
     & ,ifaceBoxhi0,ifaceBoxhi1
     & );
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifaceVello0,ifaceVello1
      integer ifaceVelhi0,ifaceVelhi1
      REAL*8 faceVel(
     & ifaceVello0:ifaceVelhi0,
     & ifaceVello1:ifaceVelhi1)
      integer icellHlo0,icellHlo1
      integer icellHhi0,icellHhi1
      REAL*8 cellH(
     & icellHlo0:icellHhi0,
     & icellHlo1:icellHhi1)
      integer ifaceHlo0,ifaceHlo1
      integer ifaceHhi0,ifaceHhi1
      REAL*8 faceH(
     & ifaceHlo0:ifaceHhi0,
     & ifaceHlo1:ifaceHhi1)
      integer ifaceDlo0,ifaceDlo1
      integer ifaceDhi0,ifaceDhi1
      REAL*8 faceD(
     & ifaceDlo0:ifaceDhi0,
     & ifaceDlo1:ifaceDhi1)
      integer dir
      REAL*8 dx
      integer ifaceBoxlo0,ifaceBoxlo1
      integer ifaceBoxhi0,ifaceBoxhi1
      integer i0,i1
      integer ii0,ii1
      REAL*8 oneOnDX
      OneOnDx = (1.0d0)/dx;
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = ifaceBoxlo1,ifaceBoxhi1
      do i0 = ifaceBoxlo0,ifaceBoxhi0
      faceVel(i0,i1) = faceVel(i0,i1)
     & + faceD(i0,i1)/faceH(i0,i1)
     & * OneOnDx
     & * ( - cellH(i0-ii0,i1-ii1) + cellH(i0,i1))
      enddo
      enddo
      return
      end
      subroutine SUBTRACTDFLUX(
     & faceFlux
     & ,ifaceFluxlo0,ifaceFluxlo1
     & ,ifaceFluxhi0,ifaceFluxhi1
     & ,cellH
     & ,icellHlo0,icellHlo1
     & ,icellHhi0,icellHhi1
     & ,faceD
     & ,ifaceDlo0,ifaceDlo1
     & ,ifaceDhi0,ifaceDhi1
     & ,dir
     & ,dx
     & ,ifaceBoxlo0,ifaceBoxlo1
     & ,ifaceBoxhi0,ifaceBoxhi1
     & );
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifaceFluxlo0,ifaceFluxlo1
      integer ifaceFluxhi0,ifaceFluxhi1
      REAL*8 faceFlux(
     & ifaceFluxlo0:ifaceFluxhi0,
     & ifaceFluxlo1:ifaceFluxhi1)
      integer icellHlo0,icellHlo1
      integer icellHhi0,icellHhi1
      REAL*8 cellH(
     & icellHlo0:icellHhi0,
     & icellHlo1:icellHhi1)
      integer ifaceDlo0,ifaceDlo1
      integer ifaceDhi0,ifaceDhi1
      REAL*8 faceD(
     & ifaceDlo0:ifaceDhi0,
     & ifaceDlo1:ifaceDhi1)
      integer dir
      REAL*8 dx
      integer ifaceBoxlo0,ifaceBoxlo1
      integer ifaceBoxhi0,ifaceBoxhi1
      integer i0,i1
      integer ii0,ii1
      REAL*8 oneOnDx
      OneOnDx = (1.0d0)/dx;
      ii0 = CHF_ID(dir,0)
                ii1 = CHF_ID(dir,1)
      do i1 = ifaceBoxlo1,ifaceBoxhi1
      do i0 = ifaceBoxlo0,ifaceBoxhi0
      faceFlux(i0,i1) = faceFlux(i0,i1)
     & + faceD(i0,i1) * OneOnDx
     & * ( - cellH(i0-ii0,i1-ii1) + cellH(i0,i1))
      enddo
      enddo
      return
      end
      subroutine VELINITIALGUESS(
     & vel
     & ,ivello0,ivello1
     & ,ivelhi0,ivelhi1
     & ,nvelcomp
     & ,rhs
     & ,irhslo0,irhslo1
     & ,irhshi0,irhshi1
     & ,nrhscomp
     & ,beta
     & ,ibetalo0,ibetalo1
     & ,ibetahi0,ibetahi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nvelcomp
      integer ivello0,ivello1
      integer ivelhi0,ivelhi1
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & 0:nvelcomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & 0:nrhscomp-1)
      integer ibetalo0,ibetalo1
      integer ibetahi0,ibetahi1
      REAL*8 beta(
     & ibetalo0:ibetahi0,
     & ibetalo1:ibetahi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
        if (beta(i0,i1) .gt. 0.001) then
           vel(i0,i1,0) = -rhs(i0,i1,0)/beta(i0,i1)
           vel(i0,i1,1) = -rhs(i0,i1,1)/beta(i0,i1)
        else
           vel(i0,i1,0) = (0.0d0)
           vel(i0,i1,1) = (0.0d0)
        endif
      enddo
      enddo
      return
      end
      subroutine UPWINDLAYERFLUX(
     & f
     & ,iflo0,iflo1
     & ,ifhi0,ifhi1
     & ,phiD
     & ,iphiDlo0,iphiDlo1
     & ,iphiDhi0,iphiDhi1
     & ,phiU
     & ,iphiUlo0,iphiUlo1
     & ,iphiUhi0,iphiUhi1
     & ,u
     & ,iulo0,iulo1
     & ,iuhi0,iuhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1)
      integer iphiDlo0,iphiDlo1
      integer iphiDhi0,iphiDhi1
      REAL*8 phiD(
     & iphiDlo0:iphiDhi0,
     & iphiDlo1:iphiDhi1)
      integer iphiUlo0,iphiUlo1
      integer iphiUhi0,iphiUhi1
      REAL*8 phiU(
     & iphiUlo0:iphiUhi0,
     & iphiUlo1:iphiUhi1)
      integer iulo0,iulo1
      integer iuhi0,iuhi1
      REAL*8 u(
     & iulo0:iuhi0,
     & iulo1:iuhi1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (u(i0,i1) . ge. (0.0d0)) then
         f(i0,i1) = u(i0,i1) * phiD(i0,i1);
      else
         f(i0,i1) = u(i0,i1) * phiU(i0,i1);
      end if
      enddo
      enddo
      return
      end
      subroutine FABMINPLUS(
     & fs
     & ,ifslo0,ifslo1
     & ,ifshi0,ifshi1
     & ,f
     & ,iflo0,iflo1
     & ,ifhi0,ifhi1
     & ,p
     & ,iplo0,iplo1
     & ,iphi0,iphi1
     & ,scale
     & ,fmax
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ifslo0,ifslo1
      integer ifshi0,ifshi1
      REAL*8 fs(
     & ifslo0:ifshi0,
     & ifslo1:ifshi1)
      integer iflo0,iflo1
      integer ifhi0,ifhi1
      REAL*8 f(
     & iflo0:ifhi0,
     & iflo1:ifhi1)
      integer iplo0,iplo1
      integer iphi0,iphi1
      REAL*8 p(
     & iplo0:iphi0,
     & iplo1:iphi1)
      REAL*8 scale
      REAL*8 fmax
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      fs(i0,i1) =
     & min(fmax, f(i0,i1) + scale * p(i0,i1))
      enddo
      enddo
      return
      end
      subroutine ADVECTFRAC_SLC_SUBOPTIMAL(
     & frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,dh
     & ,idhlo0,idhlo1
     & ,idhhi0,idhhi1
     & ,oldfrac
     & ,ioldfraclo0,ioldfraclo1
     & ,ioldfrachi0,ioldfrachi1
     & ,faceIceVel
     & ,ifaceIceVello0,ifaceIceVello1
     & ,ifaceIceVelhi0,ifaceIceVelhi1
     & ,faceCalvVel
     & ,ifaceCalvVello0,ifaceCalvVello1
     & ,ifaceCalvVelhi0,ifaceCalvVelhi1
     & ,faceIceFlux
     & ,ifaceIceFluxlo0,ifaceIceFluxlo1
     & ,ifaceIceFluxhi0,ifaceIceFluxhi1
     & ,h
     & ,ihlo0,ihlo1
     & ,ihhi0,ihhi1
     & ,dx
     & ,dt
     & ,eps
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer idhlo0,idhlo1
      integer idhhi0,idhhi1
      REAL*8 dh(
     & idhlo0:idhhi0,
     & idhlo1:idhhi1)
      integer ioldfraclo0,ioldfraclo1
      integer ioldfrachi0,ioldfrachi1
      REAL*8 oldfrac(
     & ioldfraclo0:ioldfrachi0,
     & ioldfraclo1:ioldfrachi1)
      integer ifaceIceVello0,ifaceIceVello1
      integer ifaceIceVelhi0,ifaceIceVelhi1
      REAL*8 faceIceVel(
     & ifaceIceVello0:ifaceIceVelhi0,
     & ifaceIceVello1:ifaceIceVelhi1)
      integer ifaceCalvVello0,ifaceCalvVello1
      integer ifaceCalvVelhi0,ifaceCalvVelhi1
      REAL*8 faceCalvVel(
     & ifaceCalvVello0:ifaceCalvVelhi0,
     & ifaceCalvVello1:ifaceCalvVelhi1)
      integer ifaceIceFluxlo0,ifaceIceFluxlo1
      integer ifaceIceFluxhi0,ifaceIceFluxhi1
      REAL*8 faceIceFlux(
     & ifaceIceFluxlo0:ifaceIceFluxhi0,
     & ifaceIceFluxlo1:ifaceIceFluxhi1)
      integer ihlo0,ihlo1
      integer ihhi0,ihhi1
      REAL*8 h(
     & ihlo0:ihhi0,
     & ihlo1:ihhi1)
      REAL*8 dx
      REAL*8 dt
      REAL*8 eps
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      integer i0,i1
      integer ioff0,ioff1
      REAL*8 factor, speed, onemeps, uim, ucm, uip, ucp, uc, tmp,fm,f,fp
      onemeps = 1.0d0 - eps
      factor = dt/dx
      ioff0 = CHF_ID(0,idir)
      ioff1 = CHF_ID(1,idir)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      speed = (0.0d0)
      uim = faceIceVel(i0,i1)
      ucm = faceCalvVel(i0,i1)
      uip = faceIceVel(i0+ioff0,i1+ioff1)
      ucp = faceCalvVel(i0+ioff0,i1+ioff1)
      fm = oldfrac(i0-ioff0,i1-ioff1)
      fp = oldfrac(i0+ioff0,i1+ioff1)
      f = oldfrac(i0,i1) + 1.0d-10
      if (fm .gt. onemeps) then
         speed = speed + max(uim,(0.0d0))
      end if
      if (fm .lt. eps ) then
         uc = max(ucp,max(ucm,(0.0d0)))
         speed = speed - uc
         dh(i0,i1) = dh(i0,i1)
     & - uc*factor*h(i0,i1)/f
     & - factor*min((0.0d0),faceIceFlux(i0,i1))
      end if
      if (fp .gt. onemeps) then
         speed = speed - min(uip,(0.0d0))
      end if
      if (fp .lt. eps ) then
         uc = - min(ucm,min(ucp,(0.0d0)))
         speed = speed - uc
         dh(i0,i1) = dh(i0,i1)
     & - uc*factor*h(i0,i1)/f
     & + factor*max((0.0d0),faceIceFlux(i0+ioff0,i1+ioff1))
      end if
      frac(i0,i1) = frac(i0,i1) + factor*speed
      if (frac(i0,i1) .gt. (1.0d0)) then
         frac(i0,i1) = (1.0d0)
      endif
      if (frac(i0,i1) .lt. (0.0d0)) then
         frac(i0,i1) = (0.0d0)
      endif
      enddo
      enddo
      return
      end
      subroutine ABLATERATE_SLC_SUBOPTIMAL(
     & facerate
     & ,ifaceratelo0,ifaceratelo1
     & ,ifaceratehi0,ifaceratehi1
     & ,frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,rate
     & ,iratelo0,iratelo1
     & ,iratehi0,iratehi1
     & ,eps
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifaceratelo0,ifaceratelo1
      integer ifaceratehi0,ifaceratehi1
      REAL*8 facerate(
     & ifaceratelo0:ifaceratehi0,
     & ifaceratelo1:ifaceratehi1)
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer iratelo0,iratelo1
      integer iratehi0,iratehi1
      REAL*8 rate(
     & iratelo0:iratehi0,
     & iratelo1:iratehi1)
      REAL*8 eps
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      integer i0,i1
      integer ioff0,ioff1
      integer itrans0,itrans1
      integer lu0,lu1
      integer ld0,ld1
      REAL*8 speed, grad_frac_n,grad_frac_t,fcrate
      ioff0 = CHF_ID(0,idir)
      ioff1 = CHF_ID(1,idir)
      itrans0 = CHF_ID(1,idir)
      itrans1 = CHF_ID(0,idir)
      lu0 = CHF_ID(1,idir) - CHF_ID(0,idir)
      lu1 = CHF_ID(0,idir) - CHF_ID(1,idir)
      ld0 = - CHF_ID(1,idir) - CHF_ID(0,idir)
      ld1 = - CHF_ID(0,idir) - CHF_ID(1,idir)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      speed = (0.0d0)
      fcrate = max(rate(i0,i1), rate(i0-ioff0,i1-ioff1))
      if (fcrate .gt. (0.0d0)) then
         grad_frac_t = 0.25d0 * (
     & +abs(frac(i0,i1)
     & -frac(i0+itrans0,i1+itrans1))
     & +abs(frac(i0,i1)
     & -frac(i0-itrans0,i1-itrans1))
     & +abs(frac(i0-ioff0,i1-ioff1)
     & -frac(i0+lu0,i1+lu1))
     & +abs(frac(i0-ioff0,i1-ioff1)
     & -frac(i0+ld0,i1+ld1))
     & )
         grad_frac_n = frac(i0,i1)
     & - frac(i0-ioff0,i1-ioff1)
         grad_frac_n = grad_frac_n /
     & sqrt(grad_frac_n**2 + grad_frac_t**2 + 1.0d-12)
         speed = (0.0d0)
         if (frac(i0-ioff0,i1-ioff1).lt.eps) then
            speed = speed + max(grad_frac_n,(0.0d0)) * fcrate
         end if
         if (frac(i0,i1).lt.eps) then
            speed = speed + min(grad_frac_n,(0.0d0)) * fcrate
         end if
         facerate(i0,i1) = speed
      end if
      enddo
      enddo
      return
      end
      subroutine ABLATERATECC_SLC_SUBOPTIMAL(
     & facerate
     & ,ifaceratelo0,ifaceratelo1
     & ,ifaceratehi0,ifaceratehi1
     & ,frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,rate
     & ,iratelo0,iratelo1
     & ,iratehi0,iratehi1
     & ,eps
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifaceratelo0,ifaceratelo1
      integer ifaceratehi0,ifaceratehi1
      REAL*8 facerate(
     & ifaceratelo0:ifaceratehi0,
     & ifaceratelo1:ifaceratehi1)
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer iratelo0,iratelo1
      integer iratehi0,iratehi1
      REAL*8 rate(
     & iratelo0:iratehi0,
     & iratelo1:iratehi1)
      REAL*8 eps
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      integer i0,i1
      integer ioff0,ioff1
      integer itrans0,itrans1
      REAL*8 speed, grad_frac_n,grad_frac_t
      ioff0 = CHF_ID(0,idir)
      ioff1 = CHF_ID(1,idir)
      itrans0 = CHF_ID(1,idir)
      itrans1 = CHF_ID(0,idir)
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      speed = (0.0d0)
      if ((rate(i0,i1).gt.(0.0d0))) then
         grad_frac_t = (0.500d0) * (
     & +abs(frac(i0,i1)-frac(i0+itrans0,i1+itrans1))
     & +abs(frac(i0,i1)-frac(i0-itrans0,i1-itrans1))
     & )
         speed = (0.0d0)
         if (frac(i0-ioff0,i1-ioff1).lt.eps) then
            grad_frac_n = frac(i0,i1)
     & - frac(i0-ioff0,i1-ioff1)
            grad_frac_n = grad_frac_n /
     & sqrt(grad_frac_n**2 + grad_frac_t**2 + 1.0d-12)
            speed = speed + max(grad_frac_n,(0.0d0)) * rate(i0,i1)
         end if
         if (frac(i0+ioff0,i1+ioff1).lt.eps) then
            grad_frac_n = -frac(i0,i1)
     & + frac(i0+ioff0,i1+ioff1)
            grad_frac_n = grad_frac_n /
     & sqrt(grad_frac_n**2 + grad_frac_t**2 + 1.0d-12)
            speed = speed - min(grad_frac_n,(0.0d0)) * rate(i0,i1)
         end if
         if (speed .gt. (1.0d-10)) then
            facerate(i0+ioff0,i1+ioff1) = -speed
         end if
      endif
      enddo
      enddo
      return
      end
      subroutine FINDFRONTVEL(
     & oldFrac
     & ,ioldFraclo0,ioldFraclo1
     & ,ioldFrachi0,ioldFrachi1
     & ,uFaceVel
     & ,iuFaceVello0,iuFaceVello1
     & ,iuFaceVelhi0,iuFaceVelhi1
     & ,vFaceVel
     & ,ivFaceVello0,ivFaceVello1
     & ,ivFaceVelhi0,ivFaceVelhi1
     & ,cRate
     & ,icRatelo0,icRatelo1
     & ,icRatehi0,icRatehi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,ccVel
     & ,iccVello0,iccVello1
     & ,iccVelhi0,iccVelhi1
     & ,nccVelcomp
     & ,tmpVel
     & ,itmpVello0,itmpVello1
     & ,itmpVelhi0,itmpVelhi1
     & ,ntmpVelcomp
     & ,relativeVel
     & ,irelativeVello0,irelativeVello1
     & ,irelativeVelhi0,irelativeVelhi1
     & ,nrelativeVelcomp
     & ,normalCalving
     & ,epsFrac
     & ,epsVel
     & ,dx
     & ,iStep
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ioldFraclo0,ioldFraclo1
      integer ioldFrachi0,ioldFrachi1
      REAL*8 oldFrac(
     & ioldFraclo0:ioldFrachi0,
     & ioldFraclo1:ioldFrachi1)
      integer iuFaceVello0,iuFaceVello1
      integer iuFaceVelhi0,iuFaceVelhi1
      REAL*8 uFaceVel(
     & iuFaceVello0:iuFaceVelhi0,
     & iuFaceVello1:iuFaceVelhi1)
      integer ivFaceVello0,ivFaceVello1
      integer ivFaceVelhi0,ivFaceVelhi1
      REAL*8 vFaceVel(
     & ivFaceVello0:ivFaceVelhi0,
     & ivFaceVello1:ivFaceVelhi1)
      integer icRatelo0,icRatelo1
      integer icRatehi0,icRatehi1
      REAL*8 cRate(
     & icRatelo0:icRatehi0,
     & icRatelo1:icRatehi1)
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
      integer nccVelcomp
      integer iccVello0,iccVello1
      integer iccVelhi0,iccVelhi1
      REAL*8 ccVel(
     & iccVello0:iccVelhi0,
     & iccVello1:iccVelhi1,
     & 0:nccVelcomp-1)
      integer ntmpVelcomp
      integer itmpVello0,itmpVello1
      integer itmpVelhi0,itmpVelhi1
      REAL*8 tmpVel(
     & itmpVello0:itmpVelhi0,
     & itmpVello1:itmpVelhi1,
     & 0:ntmpVelcomp-1)
      integer nrelativeVelcomp
      integer irelativeVello0,irelativeVello1
      integer irelativeVelhi0,irelativeVelhi1
      REAL*8 relativeVel(
     & irelativeVello0:irelativeVelhi0,
     & irelativeVello1:irelativeVelhi1,
     & 0:nrelativeVelcomp-1)
      integer normalCalving
      REAL*8 epsFrac
      REAL*8 epsVel
      REAL*8 dx
      integer iStep
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, ioff0,ioff1, joff0,joff1
      integer ncomp, idx(2)
      integer idir, jdir, ii, jj, isSet(2)
      integer nSheet, nIceFree, nOcean, nCorner
      integer in(2),ie(2),is(2),iw(2)
      REAL*8 vel(2), grad(2)
      REAL*8 calvingRate, critFrac
      REAL*8 u, v, speed, uUnit, vUnit
      in=0
      is=0
      ie=0
      iw=0
      critFrac=0.5d0
      ncomp = nrelativeVelcomp
      idx=0
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
        vel(1:2)=0.0d0
        isSet(1:2)=0
        idx(1:ncomp)=(/ i0,i1 /)
        nOcean=0
        nCorner=0
        if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0,i1)-1.0d0) .gt. epsFrac)) then
           do idir = 0, 2 -1
              ii=idir+1
              vel(ii)=ccVel(i0,i1,idir)
              isSet(ii)=1
           enddo
        else
           idx(1:ncomp)=(/ i0,i1 /)
           nSheet=0
           nIceFree=0
           if (abs(1.0d0-oldFrac(i0,i1)) .lt. epsFrac) then
              nSheet=1
           elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
              nIceFree=1
           endif
           do idir = 0, 2 -1
              tmpVel(i0,i1,idir)=0.0d0
              relativeVel(i0,i1,idir)=0.0d0
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
              if (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) th
     &en
                 nSheet = nSheet + 1
              elseif (abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) then
                 nIceFree = nIceFree + 1
              endif
              if (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) th
     &en
                 nSheet = nSheet + 1
              elseif (abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) then
                 nIceFree = nIceFree + 1
              endif
           enddo
           if ((nSheet .eq. 5) .or. (nIceFree .eq. 5)) then
              cycle
           endif
           if (nSheet .eq. 0) then
      ioff0=CHF_ID(0,0)
      ioff1=CHF_ID(1,0)
      joff0=CHF_ID(0,1)
      joff1=CHF_ID(1,1)
              in(1:ncomp)=(/ i0+joff0,i1+joff1 /)
              ie(1:ncomp)=(/ i0+ioff0,i1+ioff1 /)
              is(1:ncomp)=(/ i0-joff0,i1-joff1 /)
              iw(1:ncomp)=(/ i0-ioff0,i1-ioff1 /)
              if (abs(1.0d0-oldFrac(iw(1),in(2))) .lt. epsFrac) then
                 nCorner = nCorner+1
              elseif (abs(1.0d0-oldFrac(ie(1),in(2))) .lt. epsFrac) then
                 nCorner = nCorner+1
              elseif (abs(1.0d0-oldFrac(iw(1),is(2))) .lt. epsFrac) then
                 nCorner = nCorner+1
              elseif (abs(1.0d0-oldFrac(ie(1),is(2))) .lt. epsFrac) then
                 nCorner = nCorner+1
              endif
              if (nCorner .eq. 0) cycle
              if (abs(oldFrac(i0,i1)) .lt. epsFrac) then
                 if (abs(1.0d0-oldFrac(iw(1),in(2))) .lt. epsFrac) then
                    if (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. critFrac) t
     &hen
                       vel(1)=uFaceVel(i0,i1)
                       isSet(1)=1
                    elseif (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac
     &) then
                       vel(1)=ccVel(i0-ioff0,i1-ioff1,0)
                       isSet(1)=1
                    endif
                    if (abs(oldFrac(i0+joff0,i1+joff1)) .gt. critFrac) t
     &hen
                       vel(2)=vFaceVel(i0+joff0,i1+joff1)
                       isSet(2)=1
                    elseif (abs(oldFrac(i0+joff0,i1+joff1)) .gt. epsFrac
     &) then
                       vel(2)=ccVel(i0+joff0,i1+joff1,1)
                       isSet(2)=1
                    endif
                    if ((isSet(1)+isSet(2)) .eq. 1) then
                       if (isSet(1) .eq. 1) then
                          vel(2)=ccVel(i0-ioff0,i1-ioff1,1)
                          isSet(2)=1
                       else
                          vel(1)=ccVel(i0+joff0,i1+joff1,0)
                          isSet(1)=1
                       endif
                    endif
                 elseif (abs(1.0d0-oldFrac(ie(1),in(2))) .lt. epsFrac) t
     &hen
                    if (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. critFrac) t
     &hen
                       vel(1)=uFaceVel(i0+ioff0,i1+ioff1)
                       isSet(1)=1
                    elseif (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac
     &) then
                       vel(1)=ccVel(i0+ioff0,i1+ioff1,0)
                       isSet(1)=1
                    endif
                    if (abs(oldFrac(i0+joff0,i1+joff1)) .gt. critFrac) t
     &hen
                       vel(2)=vFaceVel(i0+joff0,i1+joff1)
                       isSet(2)=1
                    elseif (abs(oldFrac(i0+joff0,i1+joff1)) .gt. epsFrac
     &) then
                       vel(2)=ccVel(i0+joff0,i1+joff1,1)
                       isSet(2)=1
                    endif
                    if ((isSet(1)+isSet(2)) .eq. 1) then
                       if (isSet(1) .eq. 1) then
                          vel(2)=ccVel(i0+ioff0,i1+ioff1,1)
                          isSet(2)=1
                       else
                          vel(1)=ccVel(i0+joff0,i1+joff1,0)
                          isSet(1)=1
                       endif
                    endif
                 elseif (abs(1.0d0-oldFrac(iw(1),is(2))) .lt. epsFrac) t
     &hen
                    if (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. critFrac) t
     &hen
                       vel(1)=uFaceVel(i0,i1)
                       isSet(1)=1
                    elseif (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac
     &) then
                       vel(1)=ccVel(i0-ioff0,i1-ioff1,0)
                       isSet(1)=1
                    endif
                    if (abs(oldFrac(i0-joff0,i1-joff1)) .gt. critFrac) t
     &hen
                       vel(2)=vFaceVel(i0,i1)
                       isSet(2)=1
                    elseif (abs(oldFrac(i0-joff0,i1-joff1)) .gt. epsFrac
     &) then
                       vel(2)=ccVel(i0-joff0,i1-joff1,1)
                       isSet(2)=1
                    endif
                    if ((isSet(1)+isSet(2)) .eq. 1) then
                       if (isSet(1) .eq. 1) then
                          vel(2)=ccVel(i0-ioff0,i1-ioff1,1)
                          isSet(2)=1
                       else
                          vel(1)=ccVel(i0-joff0,i1-joff1,0)
                          isSet(1)=1
                       endif
                    endif
                 elseif (abs(1.0d0-oldFrac(ie(1),is(2))) .lt. epsFrac) t
     &hen
                    if (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. critFrac) t
     &hen
                       vel(1)=uFaceVel(i0+ioff0,i1+ioff1)
                       isSet(1)=1
                    elseif (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac
     &) then
                       vel(1)=ccVel(i0+ioff0,i1+ioff1,0)
                       isSet(1)=1
                    endif
                    if (abs(oldFrac(i0-joff0,i1-joff1)) .gt.critFrac) th
     &en
                       vel(2)=vFaceVel(i0,i1)
                       isSet(2)=1
                    elseif (abs(oldFrac(i0-joff0,i1-joff1)) .gt. epsFrac
     &) then
                       vel(2)=ccVel(i0-joff0,i1-joff1,1)
                       isSet(2)=1
                    endif
                    if ((isSet(1)+isSet(2)) .eq. 1) then
                       if (isSet(1) .eq. 1) then
                          vel(2)=ccVel(i0+ioff0,i1+ioff1,1)
                          isSet(2)=1
                       else
                          vel(1)=ccVel(i0-joff0,i1-joff1,0)
                          isSet(1)=1
                       endif
                    endif
                 endif
                 if ((isSet(1)+isSet(2)) .eq. 0) cycle
              endif
           endif
        endif
        if ((isSet(1)+isSet(2)) .eq. 0) then
           do idir = 0, 2 -1
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
              ii=idir+1
              jdir=mod(idir+1,2)
              jj=jdir+1
              if (abs(oldFrac(i0,i1)) .lt. epsFrac) then
                 if (nSheet .gt. 0) then
                    if ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsF
     &rac) .and.
     & (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac)) then
                       if (idir .eq. 0) then
                          vel(ii)=0.5d0*(uFaceVel(i0,i1) +
     & uFaceVel(i0+ioff0,i1+ioff1))
                          isSet(ii)=1
                       else
                          vel(ii)=0.5d0*(vFaceVel(i0,i1) +
     & vFaceVel(i0+ioff0,i1+ioff1))
                          isSet(ii)=1
                       endif
                    elseif (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. e
     &psFrac) then
                       if (idir .eq. 0) then
                          vel(ii)=uFaceVel(i0,i1)
                          isSet(ii)=1
                       else
                          vel(ii)=vFaceVel(i0,i1)
                          isSet(ii)=1
                       endif
                    elseif (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. e
     &psFrac) then
                       if (idir .eq. 0) then
                          vel(ii)=uFaceVel(i0+ioff0,i1+ioff1)
                          isSet(ii)=1
                       else
                          vel(ii)=vFaceVel(i0+ioff0,i1+ioff1)
                          isSet(ii)=1
                       endif
                    elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFra
     &c) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac)) then
                       vel(ii)=0.5d0*(ccVel(i0+ioff0,i1+ioff1,idir) +
     & ccVel(i0-ioff0,i1-ioff1,idir))
                       isSet(ii)=1
                    elseif (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac
     &) then
                       vel(ii)=ccVel(i0+ioff0,i1+ioff1,idir)
                       isSet(ii)=1
                    elseif (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac
     &) then
                       vel(ii)=ccVel(i0-ioff0,i1-ioff1,idir)
                       isSet(ii)=1
                    endif
                    if (isSet(jj) .eq. 0) then
                       if ((abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac
     &) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac)) then
                          vel(jj)=0.5d0*(ccVel(i0-ioff0,i1-ioff1,jdir) +
     & ccVel(i0+ioff0,i1+ioff1,jdir))
                          isSet(jj)=1
                       elseif (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsF
     &rac) then
                          vel(jj)=ccVel(i0+ioff0,i1+ioff1,jdir)
                          isSet(jj)=1
                       elseif (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsF
     &rac) then
                          vel(jj)=ccVel(i0-ioff0,i1-ioff1,jdir)
                          isSet(jj)=1
                       endif
                    endif
                 endif
              elseif (abs(1.0d0-oldFrac(i0,i1)) .lt. epsFrac) then
                 if ((nIceFree .gt. 0) .or. (nSheet .le. 3)) then
                    if ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .
     &and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac)) then
                       if (idir .eq. 0) then
                          vel(ii)=uFaceVel(i0,i1)
                          isSet(ii)=1
                       else
                          vel(ii)=vFaceVel(i0,i1)
                          isSet(ii)=1
                       endif
                    elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFra
     &c) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac)) then
                       if (idir .eq. 0) then
                          vel(ii)=uFaceVel(i0+ioff0,i1+ioff1)
                          isSet(ii)=1
                       else
                          vel(ii)=vFaceVel(i0+ioff0,i1+ioff1)
                          isSet(ii)=1
                       endif
                       isSet(ii)=1
                    else
                       vel(ii)=ccVel(i0,i1,idir)
                       isSet(ii)=1
                    endif
                 endif
              endif
           enddo
        endif
        if ((isSet(1)+isSet(2)) .eq. 0) cycle
        u=vel(1)
        v=vel(2)
        do idir = 0, 2 -1
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
           if (mask(i0+ioff0,i1+ioff1) .eq. (4)) then
              nOcean = nOcean + 1
           endif
           if (mask(i0-ioff0,i1-ioff1) .eq. (4)) then
              nOcean = nOcean + 1
           endif
        enddo
        if ((nOcean .eq. 0) .and. (abs(oldFrac(i0,i1)) .gt. epsFrac)) th
     &en
           if (mask(i0,i1) .eq. (4)) then
              nOcean = nOcean+1
           endif
      ioff0=CHF_ID(0,0)
      ioff1=CHF_ID(1,0)
      joff0=CHF_ID(0,1)
      joff1=CHF_ID(1,1)
            in(1:ncomp)=(/ i0+joff0,i1+joff1 /)
            ie(1:ncomp)=(/ i0+ioff0,i1+ioff1 /)
            is(1:ncomp)=(/ i0-joff0,i1-joff1 /)
            iw(1:ncomp)=(/ i0-ioff0,i1-ioff1 /)
            if (mask(iw(1),in(2)) .eq. (4)) then
               nOcean = nOcean+1
            elseif (mask(ie(1),in(2)) .eq. (4)) then
               nOcean = nOcean+1
            elseif (mask(iw(1),is(2)) .eq. (4)) then
               nOcean = nOcean+1
            elseif (mask(ie(1),is(2)) .eq. (4)) then
               nOcean = nOcean+1
            endif
        endif
        if (nOcean .gt. 0) then
           calvingRate=cRate(i0,i1)
           if (normalCalving .eq. 1) then
              grad(1:2)=0.0d0
              do idir = 0, 2 -1
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
                 ii=idir+1
                 if ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .and
     &.
     & (mask(i0-ioff0,i1-ioff1) .ne. (4))) then
                    grad(ii)=(oldFrac(i0+ioff0,i1+ioff1)-oldFrac(i0,i1))
     &/dx
                 elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) 
     &.and.
     & (mask(i0+ioff0,i1+ioff1) .ne. (4))) then
                    grad(ii)=(oldFrac(i0,i1)-oldFrac(i0-ioff0,i1-ioff1))
     &/dx
                 elseif ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) 
     &.and.
     & (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac)) then
                    grad(ii)=1.0d0/dx
                 elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) 
     &.and.
     & (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac)) then
                    grad(ii)=-1.0d0/dx
                 else
                    if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
                       if ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. e
     &psFrac) .and.
     & ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac))) then
                          grad(ii)=(1.0d0-0.5d0*(oldFrac(i0,i1)+oldFrac(
     &i0-ioff0,i1-ioff1)))/dx
                       elseif ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .l
     &t. epsFrac) .and.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
                          grad(ii)=(0.5d0*(oldFrac(i0,i1)+oldFrac(i0+iof
     &f0,i1+ioff1))-1.0d0)/dx
                       elseif ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. eps
     &Frac) .and.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
                          grad(ii)=0.5d0*(oldFrac(i0,i1)+oldFrac(i0+ioff
     &0,i1+ioff1))/dx
                       elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. eps
     &Frac) .and.
     & ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac))) then
                          grad(ii)=-0.5d0*(oldFrac(i0,i1)+oldFrac(i0-iof
     &f0,i1-ioff1))/dx
                       else
                          grad(ii)=0.5d0*(oldFrac(i0+ioff0,i1+ioff1)-old
     &Frac(i0-ioff0,i1-ioff1))/dx
                       endif
                    elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
                       grad(ii)=0.5d0*(oldFrac(i0+ioff0,i1+ioff1)-oldFra
     &c(i0-ioff0,i1-ioff1))/dx
                    elseif (abs(1.0d0-oldFrac(i0,i1)) .lt. epsFrac) then
                       grad(ii)=0.5d0*(oldFrac(i0+ioff0,i1+ioff1)-oldFra
     &c(i0-ioff0,i1-ioff1))/dx
                    endif
                 endif
              enddo
              speed=dsqrt(grad(1)*grad(1)+grad(2)*grad(2))
              if (speed .lt. epsVel) then
                 uUnit=0.0d0
                 vUnit=0.0d0
              else
                 uUnit=-grad(1)/speed
                 vUnit=-grad(2)/speed
              endif
           else
              speed=dsqrt(u*u+v*v)
              if (speed .lt. epsVel) then
                 uUnit=0.0d0
                 vUnit=0.0d0
              else
                 uUnit=u/speed
                 vUnit=v/speed
              endif
           endif
           tmpVel(i0,i1,0)=uUnit
           tmpVel(i0,i1,1)=vUnit
           relativeVel(i0,i1,0)=u-calvingRate*uUnit
           relativeVel(i0,i1,1)=v-calvingRate*vUnit
        else
           tmpVel(i0,i1,0)=0.0d0
           tmpVel(i0,i1,1)=0.0d0
           relativeVel(i0,i1,0)=u
           relativeVel(i0,i1,1)=v
        endif
      enddo
      enddo
      return
      end
      subroutine ADVECTFRAC(
     & frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,oldFrac
     & ,ioldFraclo0,ioldFraclo1
     & ,ioldFrachi0,ioldFrachi1
     & ,relativeFrontVel
     & ,irelativeFrontVello0,irelativeFrontVello1
     & ,irelativeFrontVelhi0,irelativeFrontVelhi1
     & ,dx
     & ,dt
     & ,epsFrac
     & ,epsVel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer ioldFraclo0,ioldFraclo1
      integer ioldFrachi0,ioldFrachi1
      REAL*8 oldFrac(
     & ioldFraclo0:ioldFrachi0,
     & ioldFraclo1:ioldFrachi1)
      integer irelativeFrontVello0,irelativeFrontVello1
      integer irelativeFrontVelhi0,irelativeFrontVelhi1
      REAL*8 relativeFrontVel(
     & irelativeFrontVello0:irelativeFrontVelhi0,
     & irelativeFrontVello1:irelativeFrontVelhi1)
      REAL*8 dx
      REAL*8 dt
      REAL*8 epsFrac
      REAL*8 epsVel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      integer i0,i1
      integer ioff0,ioff1, joff0,joff1
      integer jdir, ido
      REAL*8 factor
      REAL*8 vel
      REAL*8 lowerFrac, upperFrac
      ido=0
      factor = dt/dx
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
        jdir=mod(idir+1,2)
      joff0=CHF_ID(0,jdir)
      joff1=CHF_ID(1,jdir)
        vel=relativeFrontVel(i0,i1)
        if (abs(vel) .lt. epsVel) cycle
        lowerFrac=0.0d0
        upperFrac=0.0d0
        if ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac)) then
           lowerFrac=0.0d0
           upperFrac=0.0d0
           if (vel .gt. 0.0d0) then
              if (oldFrac(i0,i1) .gt. epsFrac) then
                 upperFrac=1.0d0
              endif
            elseif (vel .lt. 0.0d0) then
              if (oldFrac(i0,i1) .lt. (1.0d0-epsFrac)) then
                 upperFrac=1.0d0
              endif
           endif
        elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac)) then
           lowerFrac=0.0d0
           upperFrac=0.0d0
           if (vel .lt. 0.0d0) then
              if (oldFrac(i0,i1) .gt. epsFrac) then
                 lowerFrac=1.0d0
              endif
           elseif (vel .gt. 0.0d0) then
              if (oldFrac(i0,i1) .lt. (1.0d0-epsFrac)) then
                 lowerFrac=1.0d0
              endif
           endif
        elseif (((abs(oldFrac(i0-joff0,i1-joff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0+joff0,i1+joff1)) .lt. epsFrac)) .or.
     & ((abs(oldFrac(i0+joff0,i1+joff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0-joff0,i1-joff1)) .lt. epsFrac))) then
           if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
              if (vel .lt. 0.0d0) then
                 lowerFrac=oldFrac(i0,i1)
                 upperFrac=oldFrac(i0+ioff0,i1+ioff1)
              elseif (vel .gt. 0.0d0) then
                 lowerFrac=oldFrac(i0-ioff0,i1-ioff1)
                 upperFrac=oldFrac(i0,i1)
              endif
              if ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) .
     &and.
     & ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac))) then
                 lowerFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0-ioff0,i1-iof
     &f1))
                 upperFrac=1.0d0
              elseif ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFra
     &c) .and.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
                 lowerFrac=1.0d0
                 upperFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0+ioff0,i1+iof
     &f1))
              elseif ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .an
     &d.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
                 lowerFrac=0.0d0
                 upperFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0+ioff0,i1+iof
     &f1))
              elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) .an
     &d.
     & ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac))) then
                 lowerFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0-ioff0,i1-iof
     &f1))
                 upperFrac=0.0d0
              endif
           elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
              if (vel .lt. 0.0d0) then
                 if (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac)
     & then
                    lowerFrac=0.0d0
                    upperFrac=1.0d0
                 elseif (relativeFrontVel(i0+ioff0,i1+ioff1) .lt. -epsVe
     &l) then
                    lowerFrac=0.0d0
                    upperFrac=oldFrac(i0+ioff0,i1+ioff1)
                 endif
              elseif (vel .gt. 0.0d0) then
                 if (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac)
     & then
                    lowerFrac=1.0d0
                    upperFrac=0.0d0
                 elseif (relativeFrontVel(i0-ioff0,i1-ioff1) .gt. epsVel
     &) then
                    lowerFrac=oldFrac(i0-ioff0,i1-ioff1)
                    upperFrac=0.0d0
                 endif
              endif
           elseif (abs(1.0d0-oldFrac(i0,i1)) .lt. epsFrac) then
              if (vel .lt. 0.0d0) then
                 if (abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) then
                    lowerFrac=1.0d0
                    upperFrac=0.0d0
                 endif
              elseif (vel .gt. 0.0d0) then
                 if (abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) then
                    lowerFrac=0.0d0
                    upperFrac=1.0d0
                 endif
              endif
           endif
        elseif ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) .an
     &d.
     & ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac))) then
           if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
              lowerFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0-ioff0,i1-ioff1)
     &)
              upperFrac=1.0d0
           elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
              if (vel .lt. 0.0d0) then
                 lowerFrac=0.0d0
                 upperFrac=1.0d0
              endif
           endif
        elseif ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .an
     &d.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
           if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
              lowerFrac=1.0d0
              upperFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0+ioff0,i1+ioff1)
     &)
           elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
              if (vel .gt. 0.0d0) then
                 lowerFrac=1.0d0
                 upperFrac=0.0d0
              endif
           endif
        elseif ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .and.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
           if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
              lowerFrac=0.0d0
              upperFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0+ioff0,i1+ioff1)
     &)
           elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
              if ((vel .lt. -epsVel) .and.
     & (relativeFrontVel(i0+ioff0,i1+ioff1) .lt. -epsVel)) then
                 if ((abs(oldFrac(i0-joff0,i1-joff1)) .lt. epsFrac) .or.
     & (abs(oldFrac(i0+joff0,i1+joff1)) .lt. epsFrac)) then
                    lowerFrac=0.0d0
                    upperFrac=oldFrac(i0+ioff0,i1+ioff1)
                 endif
              endif
           elseif (abs(1.0d0-oldFrac(i0,i1)) .lt. epsFrac) then
              if (vel .gt. 0.0d0) then
                 lowerFrac=0.0d0
                 upperFrac=1.0d0
              endif
           endif
        elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) .and.
     & ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac))) then
           if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
              lowerFrac=0.5d0*(oldFrac(i0,i1)+oldFrac(i0-ioff0,i1-ioff1)
     &)
              upperFrac=0.0d0
           elseif (abs(oldFrac(i0,i1)) .lt. epsFrac) then
              if ((vel .gt. epsVel) .and.
     & (relativeFrontVel(i0-ioff0,i1-ioff1) .gt. epsVel)) then
                 if ((abs(oldFrac(i0-joff0,i1-joff1)) .lt. epsFrac) .or.
     & (abs(oldFrac(i0+joff0,i1+joff1)) .lt. epsFrac)) then
                    lowerFrac=oldFrac(i0-ioff0,i1-ioff1)
                    upperFrac=0.0d0
                 endif
              endif
           elseif (abs(1.0d0-oldFrac(i0,i1)) .lt. epsFrac) then
              if (vel .lt. 0.0d0) then
                 lowerFrac=1.0d0
                 upperFrac=0.0d0
              endif
           endif
        elseif (((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac) .a
     &nd.
     & (abs(oldFrac(i0-ioff0,i1-ioff1)) .gt. epsFrac)) .and.
     & ((abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .gt. epsFrac))) then
           if ((abs(oldFrac(i0,i1)) .gt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac)) then
              if (vel .lt. 0.0d0) then
                 lowerFrac=oldFrac(i0,i1)
                 upperFrac=oldFrac(i0+ioff0,i1+ioff1)
              elseif (vel .gt. 0.0d0) then
                 lowerFrac=oldFrac(i0-ioff0,i1-ioff1)
                 upperFrac=oldFrac(i0,i1)
              endif
           endif
        elseif ((abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .an
     &d.
     & (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac)) then
           if (abs(1.0d0-oldFrac(i0,i1)) .gt. epsFrac) then
              if (vel .lt. 0.0d0) then
                 lowerFrac=oldFrac(i0,i1)
                 upperFrac=1.0d0
              elseif (vel .gt. 0.0d0) then
                 lowerFrac=1.0d0
                 upperFrac=oldFrac(i0,i1)
              endif
           endif
        elseif ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .and.
     & (abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac)) then
           if (abs(oldFrac(i0,i1)) .gt. epsFrac) then
              if (vel .lt. 0.0d0) then
                 lowerFrac=oldFrac(i0,i1)
                 upperFrac=0.0d0
              elseif (vel .gt. 0.0d0) then
                 lowerFrac=0.0d0
                 upperFrac=oldFrac(i0,i1)
              endif
           endif
        endif
        frac(i0,i1) = frac(i0,i1)
     & - factor*vel*(upperFrac-lowerFrac)
        if (abs(vel*factor) .gt. 0.9d0) then
           write(*,*) "CFL violation at", i0,i1, vel*factor
        endif
      enddo
      enddo
      return
      end
      subroutine ADVECTFRACSIMPLE(
     & frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,oldFrac
     & ,ioldFraclo0,ioldFraclo1
     & ,ioldFrachi0,ioldFrachi1
     & ,relativeFrontVel
     & ,irelativeFrontVello0,irelativeFrontVello1
     & ,irelativeFrontVelhi0,irelativeFrontVelhi1
     & ,dx
     & ,dt
     & ,epsFrac
     & ,epsVel
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,idir
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer ioldFraclo0,ioldFraclo1
      integer ioldFrachi0,ioldFrachi1
      REAL*8 oldFrac(
     & ioldFraclo0:ioldFrachi0,
     & ioldFraclo1:ioldFrachi1)
      integer irelativeFrontVello0,irelativeFrontVello1
      integer irelativeFrontVelhi0,irelativeFrontVelhi1
      REAL*8 relativeFrontVel(
     & irelativeFrontVello0:irelativeFrontVelhi0,
     & irelativeFrontVello1:irelativeFrontVelhi1)
      REAL*8 dx
      REAL*8 dt
      REAL*8 epsFrac
      REAL*8 epsVel
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      integer i0,i1
      integer ioff0,ioff1, joff0,joff1
      integer jdir
      REAL*8 factor
      REAL*8 vel
      REAL*8 lowerFrac, upperFrac
      factor = dt/dx
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
        jdir=mod(idir+1,2)
      joff0=CHF_ID(0,jdir)
      joff1=CHF_ID(1,jdir)
        vel=relativeFrontVel(i0,i1)
        if (abs(vel) .lt. epsVel) cycle
        lowerFrac=0.0d0
        upperFrac=0.0d0
        if ((abs(oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac)) then
           lowerFrac=0.0d0
           upperFrac=0.0d0
           if (vel .gt. 0.0d0) then
              if (oldFrac(i0,i1) .gt. epsFrac) then
                 upperFrac=1.0d0
              endif
            elseif (vel .lt. 0.0d0) then
              if (oldFrac(i0,i1) .lt. (1.0d0-epsFrac)) then
                 upperFrac=1.0d0
              endif
           endif
        elseif ((abs(oldFrac(i0+ioff0,i1+ioff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0-ioff0,i1-ioff1)) .lt. epsFrac)) then
           lowerFrac=0.0d0
           upperFrac=0.0d0
           if (vel .lt. 0.0d0) then
              if (oldFrac(i0,i1) .gt. epsFrac) then
                 lowerFrac=1.0d0
              endif
           elseif (vel .gt. 0.0d0) then
              if (oldFrac(i0,i1) .lt. (1.0d0-epsFrac)) then
                 lowerFrac=1.0d0
              endif
           endif
        elseif ((abs(oldFrac(i0-joff0,i1-joff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0+joff0,i1+joff1)) .lt. epsFrac)) then
           if (vel .lt. 0.0d0) then
              lowerFrac=oldFrac(i0,i1)
              upperFrac=oldFrac(i0+ioff0,i1+ioff1)
           elseif (vel .gt. 0.0d0) then
              lowerFrac=oldFrac(i0-ioff0,i1-ioff1)
              upperFrac=oldFrac(i0,i1)
           endif
        elseif ((abs(oldFrac(i0+joff0,i1+joff1)) .lt. epsFrac) .and.
     & (abs(1.0d0-oldFrac(i0-joff0,i1-joff1)) .lt. epsFrac)) then
           if (vel .lt. 0.0d0) then
              lowerFrac=oldFrac(i0,i1)
              upperFrac=oldFrac(i0+ioff0,i1+ioff1)
           elseif (vel .gt. 0.0d0) then
              lowerFrac=oldFrac(i0-ioff0,i1-ioff1)
              upperFrac=oldFrac(i0,i1)
           endif
        endif
        frac(i0,i1) = frac(i0,i1)
     & - factor*vel*(upperFrac-lowerFrac)
        if (abs(vel*factor) .gt. 0.9d0) then
           write(*,*) "CFL violation at", i0,i1, vel*factor
        endif
      enddo
      enddo
      return
      end
      subroutine GETFRAC(
     & frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,oldFrac
     & ,ioldFraclo0,ioldFraclo1
     & ,ioldFrachi0,ioldFrachi1
     & ,copyFrac
     & ,icopyFraclo0,icopyFraclo1
     & ,icopyFrachi0,icopyFrachi1
     & ,relativeVel
     & ,irelativeVello0,irelativeVello1
     & ,irelativeVelhi0,irelativeVelhi1
     & ,nrelativeVelcomp
     & ,epsFrac
     & ,epsVel
     & ,iStep
     & ,icount
     & ,iOutofRange
     & ,lev
     & ,MaxIter
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer ioldFraclo0,ioldFraclo1
      integer ioldFrachi0,ioldFrachi1
      REAL*8 oldFrac(
     & ioldFraclo0:ioldFrachi0,
     & ioldFraclo1:ioldFrachi1)
      integer icopyFraclo0,icopyFraclo1
      integer icopyFrachi0,icopyFrachi1
      REAL*8 copyFrac(
     & icopyFraclo0:icopyFrachi0,
     & icopyFraclo1:icopyFrachi1)
      integer nrelativeVelcomp
      integer irelativeVello0,irelativeVello1
      integer irelativeVelhi0,irelativeVelhi1
      REAL*8 relativeVel(
     & irelativeVello0:irelativeVelhi0,
     & irelativeVello1:irelativeVelhi1,
     & 0:nrelativeVelcomp-1)
      REAL*8 epsFrac
      REAL*8 epsVel
      integer iStep
      integer icount
      integer iOutofRange
      integer lev
      integer MaxIter
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, ioff0,ioff1, joff0,joff1
      integer idir, ncomp, idx(2)
      integer in(2),ie(2),is(2),iw(2)
      integer iprint
      REAL*8 excess
      REAL*8 u, v, speed, uWeight, vWeight
      ncomp = nrelativeVelcomp
      idx = 0
      ioutOfRange=0
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
        in=0
        is=0
        ie=0
        iw=0
        idx(1:ncomp)=(/ i0,i1 /)
        iprint=0
      ioff0=CHF_ID(0,0)
      ioff1=CHF_ID(1,0)
      joff0=CHF_ID(0,1)
      joff1=CHF_ID(1,1)
        if (copyFrac(i0,i1) .gt. 1.0d0) then
           excess=(copyFrac(i0,i1)-1.0d0)
        elseif (copyFrac(i0,i1) .lt. epsFrac) then
           excess=copyFrac(i0,i1)
        else
           excess=0.0d0
        endif
        if (abs(excess) .gt. epsFrac) then
           u=relativeVel(i0,i1,0)
           v=relativeVel(i0,i1,1)
           speed=dsqrt(u*u+v*v)
           if (speed .lt. epsVel) then
              iprint=1
           endif
        endif
        frac(i0,i1) = frac(i0,i1)-excess
        if ((copyFrac(i0+ioff0,i1+ioff1) .gt. 1.0d0) .or.
     & (copyFrac(i0+ioff0,i1+ioff1) .lt. 0.0d0)) then
           u=relativeVel(i0+ioff0,i1+ioff1,0)
           if (u .lt. -epsVel) then
              v=relativeVel(i0+ioff0,i1+ioff1,1)
              uWeight=abs(u)/(abs(u)+abs(v))
              if (copyFrac(i0+ioff0,i1+ioff1) .gt. 1.0d0) then
                 excess=copyFrac(i0+ioff0,i1+ioff1)-1.0d0
              else
                 excess=copyFrac(i0+ioff0,i1+ioff1)
              endif
              frac(i0,i1)=frac(i0,i1)+excess*uWeight
           endif
        endif
        if ((copyFrac(i0-ioff0,i1-ioff1) .gt. 1.0d0) .or.
     & (copyFrac(i0-ioff0,i1-ioff1) .lt. 0.0d0)) then
           u=relativeVel(i0-ioff0,i1-ioff1,0)
           if (u .gt. epsVel) then
              v=relativeVel(i0-ioff0,i1-ioff1,1)
              uWeight=abs(u)/(abs(u)+abs(v))
              if (copyFrac(i0-ioff0,i1-ioff1) .gt. 1.0d0) then
                 excess=copyFrac(i0-ioff0,i1-ioff1)-1.0d0
              else
                 excess=copyFrac(i0-ioff0,i1-ioff1)
              endif
              frac(i0,i1)=frac(i0,i1)+excess*uWeight
           endif
        endif
        if ((copyFrac(i0+joff0,i1+joff1) .gt. 1.0d0) .or.
     & (copyFrac(i0+joff0,i1+joff1) .lt. 0.0d0)) then
           v=relativeVel(i0+joff0,i1+joff1,1)
           if (v .lt. -epsVel) then
              u=relativeVel(i0+joff0,i1+joff1,0)
              vWeight=abs(v)/(abs(u)+abs(v))
              if (copyFrac(i0+joff0,i1+joff1) .gt. 1.0d0) then
                 excess=copyFrac(i0+joff0,i1+joff1)-1.0d0
              else
                 excess=copyFrac(i0+joff0,i1+joff1)
              endif
              frac(i0,i1)=frac(i0,i1)+excess*vWeight
           endif
        endif
        if ((copyFrac(i0-joff0,i1-joff1) .gt. 1.0d0) .or.
     & (copyFrac(i0-joff0,i1-joff1) .lt. 0.0d0)) then
           v=relativeVel(i0-joff0,i1-joff1,1)
           if (v .gt. epsVel) then
              u=relativeVel(i0-joff0,i1-joff1,0)
              vWeight=abs(v)/(abs(u)+abs(v))
              if (copyFrac(i0-joff0,i1-joff1) .gt. 1.0d0) then
                 excess=copyFrac(i0-joff0,i1-joff1)-1.0d0
              else
                 excess=copyFrac(i0-joff0,i1-joff1)
              endif
              frac(i0,i1)=frac(i0,i1)+excess*vWeight
           endif
        endif
        if ((frac(i0,i1) .gt. (1.0d0+epsFrac)) .or.
     & (frac(i0,i1) .lt. -epsVel)) then
           iOutOfRange=1
           if ((icount .eq. (MaxIter-1)) .or. (iprint .eq. 1)) then
              in(1:ncomp)=(/ i0+joff0,i1+joff1 /)
              ie(1:ncomp)=(/ i0+ioff0,i1+ioff1 /)
              is(1:ncomp)=(/ i0-joff0,i1-joff1 /)
              iw(1:ncomp)=(/ i0-ioff0,i1-ioff1 /)
            if (icount .eq. (MaxIter-1)) then
              if (frac(i0,i1) .gt. (1.0d0+epsFrac)) then
                 frac(i0,i1)=1.0d0
              elseif (frac(i0,i1) .lt. -epsVel) then
                 frac(i0,i1)=0.0d0
              endif
            end if
           endif
        endif
      enddo
      enddo
      return
      end
      subroutine ISOLATEDFRAC(
     & frac
     & ,ifraclo0,ifraclo1
     & ,ifrachi0,ifrachi1
     & ,oldFrac
     & ,ioldFraclo0,ioldFraclo1
     & ,ioldFrachi0,ioldFrachi1
     & ,relativeVel
     & ,irelativeVello0,irelativeVello1
     & ,irelativeVelhi0,irelativeVelhi1
     & ,nrelativeVelcomp
     & ,iStep
     & ,lev
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifraclo0,ifraclo1
      integer ifrachi0,ifrachi1
      REAL*8 frac(
     & ifraclo0:ifrachi0,
     & ifraclo1:ifrachi1)
      integer ioldFraclo0,ioldFraclo1
      integer ioldFrachi0,ioldFrachi1
      REAL*8 oldFrac(
     & ioldFraclo0:ioldFrachi0,
     & ioldFraclo1:ioldFrachi1)
      integer nrelativeVelcomp
      integer irelativeVello0,irelativeVello1
      integer irelativeVelhi0,irelativeVelhi1
      REAL*8 relativeVel(
     & irelativeVello0:irelativeVelhi0,
     & irelativeVello1:irelativeVelhi1,
     & 0:nrelativeVelcomp-1)
      integer iStep
      integer lev
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1, ioff0,ioff1, joff0,joff1
      integer ncomp, idir
      integer nSheet(2), nOcean(2)
      integer ii, iprint
      integer in(3),ie(3),is(3),iw(3)
      REAL*8 epsFrac
      epsFrac = 1.0d-8
      ncomp = nrelativeVelcomp
      in=0
      is=0
      ie=0
      iw=0
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
        nSheet = 0
        nOcean = 0
        iprint = 0
        if ((frac(i0,i1) .lt. (1.0d0-epsFrac)) .and. (frac(i0,i1) .gt. e
     &psFrac)) then
           do idir = 0, 2 -1
      ioff0=CHF_ID(0,idir)
      ioff1=CHF_ID(1,idir)
              ii= idir+1
              if (frac(i0+ioff0,i1+ioff1) .gt. (1.0d0-epsFrac)) then
                 nSheet(ii) = nSheet(ii) + 1
              elseif (abs(frac(i0+ioff0,i1+ioff1)) .lt. epsFrac) then
                 nOcean(ii) = nOcean(ii) + 1
              endif
              if (frac(i0-ioff0,i1-ioff1) .gt. (1.0d0-epsFrac)) then
                 nSheet(ii) = nSheet(ii) + 1
              elseif (abs(frac(i0-ioff0,i1-ioff1)) .lt. epsFrac) then
                 nOcean(ii) = nOcean(ii) + 1
              endif
           enddo
           if (sum(nOcean) .eq. 0) then
              if (sum(nSheet) .ge. 3) then
                 frac(i0,i1)=1.0d0
                 iprint=1
              elseif ((nSheet(1) .eq. 2) .or. (nSheet(2) .eq. 2)) then
                 frac(i0,i1)=1.0d0
                 iprint=1
              endif
           endif
           if (sum(nSheet) .eq. 0) then
              if (sum(nOcean) .ge. 3) then
                 frac(i0,i1)=0.0d0
                 iprint=1
              elseif ((nOcean(1) .eq. 2) .or. (nOcean(2) .eq. 2)) then
                 frac(i0,i1)=0.0d0
                 iprint=1
              endif
           endif
        endif
      enddo
      enddo
      return
      end
      subroutine COMPUTEZVEL(
     & uz
     & ,iuzlo0,iuzlo1
     & ,iuzhi0,iuzhi1
     & ,nuzcomp
     & ,uzs
     & ,iuzslo0,iuzslo1
     & ,iuzshi0,iuzshi1
     & ,ux
     & ,iuxlo0,iuxlo1
     & ,iuxhi0,iuxhi1
     & ,nuxcomp
     & ,uy
     & ,iuylo0,iuylo1
     & ,iuyhi0,iuyhi1
     & ,nuycomp
     & ,divuhxy
     & ,idivuhxylo0,idivuhxylo1
     & ,idivuhxyhi0,idivuhxyhi1
     & ,ndivuhxycomp
     & ,fsig
     & ,ifsighi0
     & ,csig
     & ,icsighi0
     & ,dsig
     & ,idsighi0
     & ,dsx
     & ,idsxlo0,idsxlo1
     & ,idsxhi0,idsxhi1
     & ,dhx
     & ,idhxlo0,idhxlo1
     & ,idhxhi0,idhxhi1
     & ,dsy
     & ,idsylo0,idsylo1
     & ,idsyhi0,idsyhi1
     & ,dhy
     & ,idhylo0,idhylo1
     & ,idhyhi0,idhyhi1
     & ,dst
     & ,idstlo0,idstlo1
     & ,idsthi0,idsthi1
     & ,dht
     & ,idhtlo0,idhtlo1
     & ,idhthi0,idhthi1
     & ,smb
     & ,ismblo0,ismblo1
     & ,ismbhi0,ismbhi1
     & ,bmb
     & ,ibmblo0,ibmblo1
     & ,ibmbhi0,ibmbhi1
     & ,nlay
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nuzcomp
      integer iuzlo0,iuzlo1
      integer iuzhi0,iuzhi1
      REAL*8 uz(
     & iuzlo0:iuzhi0,
     & iuzlo1:iuzhi1,
     & 0:nuzcomp-1)
      integer iuzslo0,iuzslo1
      integer iuzshi0,iuzshi1
      REAL*8 uzs(
     & iuzslo0:iuzshi0,
     & iuzslo1:iuzshi1)
      integer nuxcomp
      integer iuxlo0,iuxlo1
      integer iuxhi0,iuxhi1
      REAL*8 ux(
     & iuxlo0:iuxhi0,
     & iuxlo1:iuxhi1,
     & 0:nuxcomp-1)
      integer nuycomp
      integer iuylo0,iuylo1
      integer iuyhi0,iuyhi1
      REAL*8 uy(
     & iuylo0:iuyhi0,
     & iuylo1:iuyhi1,
     & 0:nuycomp-1)
      integer ndivuhxycomp
      integer idivuhxylo0,idivuhxylo1
      integer idivuhxyhi0,idivuhxyhi1
      REAL*8 divuhxy(
     & idivuhxylo0:idivuhxyhi0,
     & idivuhxylo1:idivuhxyhi1,
     & 0:ndivuhxycomp-1)
      integer ifsighi0
      REAL*8 fsig(
     & 0:ifsighi0)
      integer icsighi0
      REAL*8 csig(
     & 0:icsighi0)
      integer idsighi0
      REAL*8 dsig(
     & 0:idsighi0)
      integer idsxlo0,idsxlo1
      integer idsxhi0,idsxhi1
      REAL*8 dsx(
     & idsxlo0:idsxhi0,
     & idsxlo1:idsxhi1)
      integer idhxlo0,idhxlo1
      integer idhxhi0,idhxhi1
      REAL*8 dhx(
     & idhxlo0:idhxhi0,
     & idhxlo1:idhxhi1)
      integer idsylo0,idsylo1
      integer idsyhi0,idsyhi1
      REAL*8 dsy(
     & idsylo0:idsyhi0,
     & idsylo1:idsyhi1)
      integer idhylo0,idhylo1
      integer idhyhi0,idhyhi1
      REAL*8 dhy(
     & idhylo0:idhyhi0,
     & idhylo1:idhyhi1)
      integer idstlo0,idstlo1
      integer idsthi0,idsthi1
      REAL*8 dst(
     & idstlo0:idsthi0,
     & idstlo1:idsthi1)
      integer idhtlo0,idhtlo1
      integer idhthi0,idhthi1
      REAL*8 dht(
     & idhtlo0:idhthi0,
     & idhtlo1:idhthi1)
      integer ismblo0,ismblo1
      integer ismbhi0,ismbhi1
      REAL*8 smb(
     & ismblo0:ismbhi0,
     & ismblo1:ismbhi1)
      integer ibmblo0,ibmblo1
      integer ibmbhi0,ibmbhi1
      REAL*8 bmb(
     & ibmblo0:ibmbhi0,
     & ibmblo1:ibmbhi1)
      integer nlay
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 coldivuhxy(0:nlay-1)
      REAL*8 colux(0:nlay), coluy(0:nlay)
      REAL*8 coluz(0:nlay)
      integer layer
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      do layer = 0, nlay - 1
         coldivuhxy(layer) = divuhxy(i0,i1, layer)
      end do
      do layer = 0, nlay
         colux(layer) = ux(i0,i1, layer)
         coluy(layer) = uy(i0,i1, layer)
      end do
      call column_compute_z_vel(coluz, uzs(i0,i1),
     & colux, coluy, coldivuhxy, fsig, dsig,
     & dsx(i0,i1), dhx(i0,i1),
     & dsy(i0,i1), dhy(i0,i1),
     & dst(i0,i1), dht(i0,i1),
     & smb(i0,i1), bmb(i0,i1))
      do layer = 0, nlay
         uz(i0,i1, layer) = coluz(layer)
      end do
      enddo
      enddo
      return
      end
      subroutine COMPUTESIGMAVEL(
     & usig
     & ,iusiglo0,iusiglo1
     & ,iusighi0,iusighi1
     & ,nusigcomp
     & ,ux
     & ,iuxlo0,iuxlo1
     & ,iuxhi0,iuxhi1
     & ,nuxcomp
     & ,uy
     & ,iuylo0,iuylo1
     & ,iuyhi0,iuyhi1
     & ,nuycomp
     & ,divuhxy
     & ,idivuhxylo0,idivuhxylo1
     & ,idivuhxyhi0,idivuhxyhi1
     & ,ndivuhxycomp
     & ,dsig
     & ,idsighi0
     & ,dht
     & ,idhtlo0,idhtlo1
     & ,idhthi0,idhthi1
     & ,smb
     & ,ismblo0,ismblo1
     & ,ismbhi0,ismbhi1
     & ,bmb
     & ,ibmblo0,ibmblo1
     & ,ibmbhi0,ibmbhi1
     & ,nlay
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nusigcomp
      integer iusiglo0,iusiglo1
      integer iusighi0,iusighi1
      REAL*8 usig(
     & iusiglo0:iusighi0,
     & iusiglo1:iusighi1,
     & 0:nusigcomp-1)
      integer nuxcomp
      integer iuxlo0,iuxlo1
      integer iuxhi0,iuxhi1
      REAL*8 ux(
     & iuxlo0:iuxhi0,
     & iuxlo1:iuxhi1,
     & 0:nuxcomp-1)
      integer nuycomp
      integer iuylo0,iuylo1
      integer iuyhi0,iuyhi1
      REAL*8 uy(
     & iuylo0:iuyhi0,
     & iuylo1:iuyhi1,
     & 0:nuycomp-1)
      integer ndivuhxycomp
      integer idivuhxylo0,idivuhxylo1
      integer idivuhxyhi0,idivuhxyhi1
      REAL*8 divuhxy(
     & idivuhxylo0:idivuhxyhi0,
     & idivuhxylo1:idivuhxyhi1,
     & 0:ndivuhxycomp-1)
      integer idsighi0
      REAL*8 dsig(
     & 0:idsighi0)
      integer idhtlo0,idhtlo1
      integer idhthi0,idhthi1
      REAL*8 dht(
     & idhtlo0:idhthi0,
     & idhtlo1:idhthi1)
      integer ismblo0,ismblo1
      integer ismbhi0,ismbhi1
      REAL*8 smb(
     & ismblo0:ismbhi0,
     & ismblo1:ismbhi1)
      integer ibmblo0,ibmblo1
      integer ibmbhi0,ibmbhi1
      REAL*8 bmb(
     & ibmblo0:ibmbhi0,
     & ibmblo1:ibmbhi1)
      integer nlay
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      REAL*8 coldivuhxy(0:nlay-1)
      REAL*8 colux(0:nlay), coluy(0:nlay)
      REAL*8 colusig(0:nlay)
      integer layer
      integer i0,i1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      do layer = 0, nlay - 1
         coldivuhxy(layer) = divuhxy(i0,i1, layer)
      end do
      do layer = 0, nlay
         colux(layer) = ux(i0,i1, layer)
         coluy(layer) = uy(i0,i1, layer)
      end do
      call column_compute_sigma_vel( colusig, colux, coluy, coldivuhxy,
     & dsig, nlay, dht(i0,i1),
     & smb(i0,i1), bmb(i0,i1))
      do layer = 0, nlay
         usig(i0,i1, layer) = colusig(layer)
      end do
      enddo
      enddo
      return
      end
      subroutine EVOLVEGROUNDEDBED(
     & newh
     & ,inewhlo0,inewhlo1
     & ,inewhhi0,inewhhi1
     & ,oldh
     & ,ioldhlo0,ioldhlo1
     & ,ioldhhi0,ioldhhi1
     & ,topg
     & ,itopglo0,itopglo1
     & ,itopghi0,itopghi1
     & ,mask
     & ,imasklo0,imasklo1
     & ,imaskhi0,imaskhi1
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer inewhlo0,inewhlo1
      integer inewhhi0,inewhhi1
      REAL*8 newh(
     & inewhlo0:inewhhi0,
     & inewhlo1:inewhhi1)
      integer ioldhlo0,ioldhlo1
      integer ioldhhi0,ioldhhi1
      REAL*8 oldh(
     & ioldhlo0:ioldhhi0,
     & ioldhlo1:ioldhhi1)
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
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      REAL*8 dh
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
      if (mask (i0,i1) .eq. (1)) then
        dh = newh(i0,i1) - oldh(i0,i1)
        topg(i0,i1) = topg(i0,i1) - dh
      end if
      enddo
      enddo
      return
      end
