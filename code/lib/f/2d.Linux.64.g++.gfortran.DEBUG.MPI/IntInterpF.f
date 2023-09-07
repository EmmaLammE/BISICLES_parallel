      subroutine INTINTERPCONSTANT(
     & fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,iblo0,iblo1
     & ,ibhi0,ibhi1
     & ,ref_ratio
     & ,ibreflo0,ibreflo1
     & ,ibrefhi0,ibrefhi1
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      integer fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      integer coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer ref_ratio
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer var
      integer ic0,ic1
      integer if0,if1
      integer ii0,ii1
      do var = 0, ncoarsecomp - 1
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0
               if0 = ic0*ref_ratio + ii0
               if1 = ic1*ref_ratio + ii1
               fine(if0,if1,var) = coarse(ic0,ic1,var)
      enddo
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine PROLONGQUAD_ICE(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,ifineBoxlo0,ifineBoxlo1
     & ,ifineBoxhi0,ifineBoxhi1
     & ,refRatio
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer ifineBoxlo0,ifineBoxlo1
      integer ifineBoxhi0,ifineBoxhi1
      integer refRatio
      INTEGER ncomp, n
      integer i ,j
      integer ic,jc
      REAL*8 L12,R12,L14,R14,L24,R24,C14,C24
      data L12,R12/0.125d0,-0.125d0/
      data L14,R14/0.21875d0,-0.15625d0/
      data L24,R24/0.03125d0,-0.09375d0/
      data C14,C24 /-0.0625d0,0.0625d0/
      ncomp = nphicomp
      do n = 0, ncomp-1
      do j = ifineBoxlo1,ifineBoxhi1
      do i = ifineBoxlo0,ifineBoxhi0
         ic = (i+1048576)/refRatio - 1048576/refRatio
         jc = (j+1048576)/refRatio - 1048576/refRatio
         phi(i,j,n) = phi(i,j,n) +
     & coarse(ic,jc,n)
         if( refRatio == 2 ) then
            if ( ic*2.lt.i ) then
               phi(i,j,n) = phi(i,j,n) +
     $ L12*coarse(ic+1,jc,n) +
     $ R12*coarse(ic-1,jc,n)
            else
               phi(i,j,n) = phi(i,j,n) +
     $ L12*coarse(ic-1,jc,n) +
     $ R12*coarse(ic+1,jc,n)
            endif
            if ( jc*2.lt.j ) then
               phi(i,j,n) = phi(i,j,n) +
     $ L12*coarse(ic,jc+1,n) +
     $ R12*coarse(ic,jc-1,n)
            else
               phi(i,j,n) = phi(i,j,n) +
     $ L12*coarse(ic,jc-1,n) +
     $ R12*coarse(ic,jc+1,n)
            endif
          else if( refRatio == 4 ) then
             if ( i - ic*4 .lt. 4/2) then
                if( (i-ic*4) == 0 ) then
                   phi(i,j,n) = phi(i,j,n)
     $ + L14*coarse(ic-1,jc,n)
     $ + (C14-1.d0)*coarse(ic,jc,n)
     $ + R14*coarse(ic+1,jc,n)
                else
                   phi(i,j,n) = phi(i,j,n)
     $ + L24*coarse(ic-1,jc,n)
     $ + (C24-1.d0)*coarse(ic,jc,n)
     $ + R24*coarse(ic+1,jc,n)
                endif
             else
                if( (i-ic*4) == 2 ) then
                   phi(i,j,n) = phi(i,j,n)
     $ + L24*coarse(ic+1,jc,n)
     $ + (C24-1.d0)*coarse(ic,jc,n)
     $ + R24*coarse(ic-1,jc,n)
                else
                   phi(i,j,n) = phi(i,j,n)
     $ + L14*coarse(ic+1,jc,n)
     $ + (C14-1.d0)*coarse(ic,jc,n)
     $ + R14*coarse(ic-1,jc,n)
                endif
             endif
             if ( j - jc*4 .lt. 4/2 ) then
                if( (j-jc*4) == 0 ) then
                   phi(i,j,n) = phi(i,j,n)
     $ + L14*coarse(ic,jc-1,n)
     $ + (C14-1.d0)*coarse(ic,jc,n)
     $ + R14*coarse(ic,jc+1,n)
                else
                   phi(i,j,n) = phi(i,j,n)
     $ + L24*coarse(ic,jc-1,n)
     $ + (C24-1.d0)*coarse(ic,jc,n)
     $ + R24*coarse(ic,jc+1,n)
                endif
             else
                if( (j-jc*4) == 2 ) then
                   phi(i,j,n) = phi(i,j,n)
     $ + L24*coarse(ic,jc+1,n)
     $ + (C24-1.d0)*coarse(ic,jc,n)
     $ + R24*coarse(ic,jc-1,n)
                else
                   phi(i,j,n) = phi(i,j,n)
     $ + L14*coarse(ic,jc+1,n)
     $ + (C14-1.d0)*coarse(ic,jc,n)
     $ + R14*coarse(ic,jc-1,n)
                endif
             endif
          endif
      enddo
      enddo
       enddo
       return
       end
      subroutine PROLONGLINEAR_ICE(
     & phi
     & ,iphilo0,iphilo1
     & ,iphihi0,iphihi1
     & ,nphicomp
     & ,coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,ifineBoxlo0,ifineBoxlo1
     & ,ifineBoxhi0,ifineBoxhi1
     & ,refRatio
     & )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & 0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer ifineBoxlo0,ifineBoxlo1
      integer ifineBoxhi0,ifineBoxhi1
      integer refRatio
      INTEGER ncomp, n
      integer i ,j
      integer ic,jc
      ncomp = nphicomp
      do n = 0, ncomp-1
      do j = ifineBoxlo1,ifineBoxhi1
      do i = ifineBoxlo0,ifineBoxhi0
         ic = (i+1048576)/refRatio - 1048576/refRatio
         jc = (j+1048576)/refRatio - 1048576/refRatio
         phi(i,j,n) = phi(i,j,n) +
     & coarse(ic,jc,n)
         if ( ic*refRatio.lt.i ) then
            phi(i,j,n) = phi(i,j,n) +
     & (coarse(ic+1,jc,n)
     & - coarse(ic,jc,n))/refRatio*(i+(0.500d0)-ic*refRatio-(0.500d0)*re
     &fRatio)
         else
            phi(i,j,n) = phi(i,j,n) +
     & (- coarse(ic-1,jc,n)
     & + coarse(ic,jc,n))/refRatio*(i+(0.500d0)-ic*refRatio-(0.500d0)*re
     &fRatio)
         endif
         if ( jc*refRatio.lt.j ) then
            phi(i,j,n) = phi(i,j,n) +
     & (coarse(ic,jc+1,n)
     & - coarse(ic,jc,n))/refRatio*(j+(0.500d0)-jc*refRatio-(0.500d0)*re
     &fRatio)
         else
            phi(i,j,n) = phi(i,j,n) +
     & (- coarse(ic,jc-1,n)
     & + coarse(ic,jc,n))/refRatio*(j+(0.500d0)-jc*refRatio-(0.500d0)*re
     &fRatio)
         endif
      enddo
      enddo
       enddo
       return
       end
