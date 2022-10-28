      subroutine GSRBHELMHOLTZAVC2D(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0,iaCoeflo1
     &           ,iaCoefhi0,iaCoefhi1
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0,ibCoef0lo1
     &           ,ibCoef0hi0,ibCoef0hi1
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0,ibCoef1lo1
     &           ,ibCoef1hi0,ibCoef1hi1
     &           ,nbCoef1comp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1
     &           ,ilambdahi0,ilambdahi1
     &           ,icol
     &           ,jcol
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 alpha
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           iaCoeflo1:iaCoefhi1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           ibCoef0lo1:ibCoef0hi1,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           ibCoef1lo1:ibCoef1hi1,
     &           0:nbCoef1comp-1)
      integer ilambdalo0,ilambdalo1
      integer ilambdahi0,ilambdahi1
      REAL*8 lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1)
      integer icol
      integer jcol
      REAL*8 dxinv,lofphi,quarter
      integer i,j
      if (2 .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      if (2 .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = (1.0d0)/(dx*dx)
      quarter = (0.500d0) * (0.500d0) *  dxinv
         do j=iregionlo1 + jcol, iregionhi1,3
            do i = iregionlo0 + icol, iregionhi0, 3
              lofphi =
     &            alpha * aCoef(i,j) * phi(i,j)
     &          - beta  *
     &             (
     &               bCoef0(i+1,j  ,0)
     &               * (phi(i+1,j  )-phi(i  ,j  ))
     &             - bCoef0(i  ,j  ,0)
     &               * (phi(i  ,j  )-phi(i-1,j  )) 
     &             + bCoef1(i  ,j+1,1)
     &               * (phi(i  ,j+1)-phi(i  ,j  ))
     &             - bCoef1(i  ,j  ,1)
     &               * (phi(i  ,j  )-phi(i  ,j-1)) 
     &             ) * dxinv
              lofphi = lofphi - beta * (
     &             + quarter * bCoef0(i+1,j  ,1)
     &             * (phi(i+1,j+1  ) -phi(i+1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             - quarter * bCoef0(i,j  ,1)
     &             * (phi(i-1,j+1  ) -phi(i-1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             + quarter * bCoef1(i  ,j+1  ,0)
     &               * (phi(i+1  ,j+1  )-phi(i-1,j+1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) )
     &             - quarter * bCoef1(i  ,j  ,0)
     &               * (phi(i+1  ,j-1  )-phi(i-1,j-1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) ) )            
              phi(i,j) = phi(i,j)
     &          - lambda(i,j) * (lofphi - rhs(i,j))
           enddo
         enddo
      return
      end
      subroutine AVCCOMPUTEOP2D(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0,iaCoeflo1
     &           ,iaCoefhi0,iaCoefhi1
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0,ibCoef0lo1
     &           ,ibCoef0hi0,ibCoef0hi1
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0,ibCoef1lo1
     &           ,ibCoef1hi0,ibCoef1hi1
     &           ,nbCoef1comp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL*8 alpha
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           iaCoeflo1:iaCoefhi1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           ibCoef0lo1:ibCoef0hi1,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           ibCoef1lo1:ibCoef1hi1,
     &           0:nbCoef1comp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dxinv, quarter
      integer i,j
      if (2 .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      if (2 .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = (1.0d0)/(dx*dx)
      quarter = (0.500d0) * (0.500d0) * dxinv
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          lofphi(i,j) =
     &       alpha * aCoef(i,j) * phi(i,j)
     &       - beta  *
     &         (
     &           bCoef0(i+1,j  ,0)
     &           * (phi(i+1,j  ) - phi(i  ,j  ))
     &         - bCoef0(i  ,j  ,0)
     &           * (phi(i  ,j  ) - phi(i-1,j  ))
     &         + bCoef1(i  ,j+1,1)
     &           * (phi(i  ,j+1) - phi(i  ,j  ))
     &         - bCoef1(i  ,j  ,1)
     &           * (phi(i  ,j  ) - phi(i  ,j-1)) 
     &         ) * dxinv
      lofphi(i,j) = lofphi(i,j) - beta * (
     &             + quarter * bCoef0(i+1,j  ,1)
     &             * (phi(i+1,j+1  ) -phi(i+1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             - quarter * bCoef0(i,j  ,1)
     &             * (phi(i-1,j+1  ) -phi(i-1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             + quarter * bCoef1(i  ,j+1  ,0)
     &               * (phi(i+1  ,j+1  )-phi(i-1,j+1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) )
     &             - quarter * bCoef1(i  ,j  ,0)
     &               * (phi(i+1  ,j-1  )-phi(i-1,j-1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) ))            
      enddo
      enddo
      return
      end
      subroutine AVCCOMPUTERES2D(
     &           res
     &           ,ireslo0,ireslo1
     &           ,ireshi0,ireshi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0,iaCoeflo1
     &           ,iaCoefhi0,iaCoefhi1
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0,ibCoef0lo1
     &           ,ibCoef0hi0,ibCoef0hi1
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0,ibCoef1lo1
     &           ,ibCoef1hi0,ibCoef1hi1
     &           ,nbCoef1comp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL*8 res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1)
      REAL*8 alpha
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           iaCoeflo1:iaCoefhi1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           ibCoef0lo1:ibCoef0hi1,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           ibCoef1lo1:ibCoef1hi1,
     &           0:nbCoef1comp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 dxinv, quarter
      integer i,j
      if (2 .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      if (2 .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = (1.0d0)/(dx*dx)
      quarter = (0.500d0) * (0.500d0) * dxinv
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          res(i,j) =
     &        rhs(i,j)
     &      - (alpha * aCoef(i,j) * phi(i,j)
     &       - beta  *
     &          (
     &            bCoef0(i+1,j  ,0)
     &            * (phi(i+1,j  ) - phi(i  ,j  ))
     &          - bCoef0(i  ,j  ,0)
     &            * (phi(i  ,j  ) - phi(i-1,j  )) 
     &          + bCoef1(i  ,j+1,1)
     &            * (phi(i  ,j+1) - phi(i  ,j  ))
     &          - bCoef1(i  ,j  ,1)
     &            * (phi(i  ,j  ) - phi(i  ,j-1)) 
     &          ) * dxinv
     &        )
          res(i,j) = res(i,j) + beta * (
     &             + quarter * bCoef0(i+1,j  ,1)
     &             * (phi(i+1,j+1  ) -phi(i+1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             - quarter * bCoef0(i,j  ,1)
     &             * (phi(i-1,j+1  ) -phi(i-1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             + quarter * bCoef1(i  ,j+1  ,0)
     &               * (phi(i+1  ,j+1  )-phi(i-1,j+1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) )
     &             - quarter * bCoef1(i  ,j  ,0)
     &               * (phi(i+1  ,j-1  )-phi(i-1,j-1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) ))            
      enddo
      enddo
      return
      end
      subroutine RESTRICTRESAVC2D(
     &           res
     &           ,ireslo0,ireslo1
     &           ,ireshi0,ireshi1
     &           ,phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0,iaCoeflo1
     &           ,iaCoefhi0,iaCoefhi1
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0,ibCoef0lo1
     &           ,ibCoef0hi0,ibCoef0hi1
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0,ibCoef1lo1
     &           ,ibCoef1hi0,ibCoef1hi1
     &           ,nbCoef1comp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,dx
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ireslo0,ireslo1
      integer ireshi0,ireshi1
      REAL*8 res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1)
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1)
      REAL*8 alpha
      integer iaCoeflo0,iaCoeflo1
      integer iaCoefhi0,iaCoefhi1
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           iaCoeflo1:iaCoefhi1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1
      integer ibCoef0hi0,ibCoef0hi1
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           ibCoef0lo1:ibCoef0hi1,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1
      integer ibCoef1hi0,ibCoef1hi1
      REAL*8 bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           ibCoef1lo1:ibCoef1hi1,
     &           0:nbCoef1comp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi,quarter
      integer i,j
      integer ii,jj
      dxinv = (1.0d0) / (dx*dx)
      denom = 2  *2
      quarter = (0.500d0) * (0.500d0) * dxinv
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          ii = i/2 
          jj = j/2 
          lofphi =
     &        alpha * aCoef(i,j) * phi(i,j)
     &      - beta  *
     &         (
     &           bCoef0(i+1,j  ,0)
     &           * (phi(i+1,j  )-phi(i  ,j  ))
     &         - bCoef0(i  ,j  ,0)
     &           * (phi(i  ,j  )-phi(i-1,j  )) 
     &         + bCoef1(i  ,j+1,1)
     &           * (phi(i  ,j+1)-phi(i  ,j  ))
     &         - bCoef1(i  ,j  ,1)
     &           * (phi(i  ,j  )-phi(i  ,j-1)) 
     &         ) * dxinv
              lofphi = lofphi - beta * (
     &             + quarter * bCoef0(i+1,j  ,1)
     &             * (phi(i+1,j+1  ) -phi(i+1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             - quarter * bCoef0(i,j  ,1)
     &             * (phi(i-1,j+1  ) -phi(i-1  ,j-1  )
     &             + phi(i,j+1  ) - phi(i,j-1  ) )
     &             + quarter * bCoef1(i  ,j+1  ,0)
     &               * (phi(i+1  ,j+1  )-phi(i-1,j+1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) )
     &             - quarter * bCoef1(i  ,j  ,0)
     &               * (phi(i+1  ,j-1  )-phi(i-1,j-1  ) 
     &             + phi(i+1  ,j  )-phi(i-1,j  ) ) )            
          res(ii,jj) = res(ii,jj)
     &                            + (rhs(i,j) - lofphi) / denom
      enddo
      enddo
      return
      end
      subroutine SUMAFACES(
     &           lhs
     &           ,ilhslo0,ilhslo1
     &           ,ilhshi0,ilhshi1
     &           ,beta
     &           ,bCoefs
     &           ,ibCoefslo0,ibCoefslo1
     &           ,ibCoefshi0,ibCoefshi1
     &           ,nbCoefscomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,dir
     &           ,scale
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ilhslo0,ilhslo1
      integer ilhshi0,ilhshi1
      REAL*8 lhs(
     &           ilhslo0:ilhshi0,
     &           ilhslo1:ilhshi1)
      REAL*8 beta
      integer nbCoefscomp
      integer ibCoefslo0,ibCoefslo1
      integer ibCoefshi0,ibCoefshi1
      REAL*8 bCoefs(
     &           ibCoefslo0:ibCoefshi0,
     &           ibCoefslo1:ibCoefshi1,
     &           0:nbCoefscomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer dir
      REAL*8 scale
      REAL*8 sumVal,quarter
      integer i,j
      integer ii,jj
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
      sumVal = bCoefs(i+ii,j+jj,dir)
     &     + bCoefs(i   ,j   ,dir)
      lhs(i,j) = lhs(i,j) + scale * beta * sumVal
      enddo
      enddo
      return
      end
