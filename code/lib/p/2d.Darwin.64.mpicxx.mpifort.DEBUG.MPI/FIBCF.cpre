      subroutine FILLHOLES(
     &           data
     &           ,idatalo0,idatalo1
     &           ,idatahi0,idatahi1
     &           ,holeVal
     &           ,eps
     &           ,numNeigbor
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idatalo0,idatalo1
      integer idatahi0,idatahi1
      REAL*8 data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1)
      REAL*8 holeVal
      REAL*8 eps
      integer numNeigbor
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ii0,ii1
      integer i0,i1
      integer count
      REAL*8 sum
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
         if (abs(data(i0,i1) - holeVal).lt.eps) then
            write(6,*) 'in conditional'
         endif
      enddo
      enddo
      return
      end
      subroutine NODETOCELL(
     &           nodeData
     &           ,inodeDatalo0,inodeDatalo1
     &           ,inodeDatahi0,inodeDatahi1
     &           ,cellData
     &           ,icellDatalo0,icellDatalo1
     &           ,icellDatahi0,icellDatahi1
     &           ,icellBoxlo0,icellBoxlo1
     &           ,icellBoxhi0,icellBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer inodeDatalo0,inodeDatalo1
      integer inodeDatahi0,inodeDatahi1
      REAL*8 nodeData(
     &           inodeDatalo0:inodeDatahi0,
     &           inodeDatalo1:inodeDatahi1)
      integer icellDatalo0,icellDatalo1
      integer icellDatahi0,icellDatahi1
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0,
     &           icellDatalo1:icellDatahi1)
      integer icellBoxlo0,icellBoxlo1
      integer icellBoxhi0,icellBoxhi1
      integer i,j
      REAL*8 w,t
      w=(0.500d0)
      w=w*(0.500d0)
      do j = icellBoxlo1,icellBoxhi1
      do i = icellBoxlo0,icellBoxhi0
      t = nodeData(i,j)
     &     + nodeData(i+1,j)
      t = t + nodeData(i,j+1)
     &     + nodeData(i+1,j+1) 
      cellData(i,j) = w * t
      enddo
      enddo
      return 
      end
      subroutine NODETOCELLCISMVEL(
     &           nodeData
     &           ,inodeDatalo0,inodeDatalo1
     &           ,inodeDatahi0,inodeDatahi1
     &           ,nnodeDatacomp
     &           ,cellData
     &           ,icellDatalo0,icellDatalo1
     &           ,icellDatahi0,icellDatahi1
     &           ,ncellDatacomp
     &           ,nGhost
     &           ,offset
     &           ,icellBoxlo0,icellBoxlo1
     &           ,icellBoxhi0,icellBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nnodeDatacomp
      integer inodeDatalo0,inodeDatalo1
      integer inodeDatahi0,inodeDatahi1
      REAL*8 nodeData(
     &           inodeDatalo0:inodeDatahi0,
     &           inodeDatalo1:inodeDatahi1,
     &           0:nnodeDatacomp-1)
      integer ncellDatacomp
      integer icellDatalo0,icellDatalo1
      integer icellDatahi0,icellDatahi1
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0,
     &           icellDatalo1:icellDatahi1,
     &           0:ncellDatacomp-1)
      integer nGhost
      integer offset
      integer icellBoxlo0,icellBoxlo1
      integer icellBoxhi0,icellBoxhi1
      integer i,j
      integer nlayers, n
      integer verticalOffset
      REAL*8 w,t
      w=(0.500d0)
      w=w*(0.500d0)
      nlayers = (icellDatahi0 - icellDatalo0 +1)
      verticalOffset = -icellDatalo0
      verticalOffset = -icellDatalo1
      do n=0,nlayers-1
      do j = icellBoxlo1,icellBoxhi1
      do i = icellBoxlo0,icellBoxhi0
      t = nodeData(n,i,j+verticalOffset)
     &     + nodeData(n,i+1,j+verticalOffset)
      t = t + nodeData(n,i,j+1+verticalOffset)
     &     + nodeData(n,i+1,j+1+verticalOffset) 
      cellData(n,i,j+verticalOffset) = w * t
      enddo
      enddo
      enddo
      return 
      end
      subroutine CELLTONODE(
     &           nodeData
     &           ,inodeDatalo0,inodeDatalo1
     &           ,inodeDatahi0,inodeDatahi1
     &           ,cellData
     &           ,icellDatalo0,icellDatalo1
     &           ,icellDatahi0,icellDatahi1
     &           ,inodeBoxlo0,inodeBoxlo1
     &           ,inodeBoxhi0,inodeBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer inodeDatalo0,inodeDatalo1
      integer inodeDatahi0,inodeDatahi1
      REAL*8 nodeData(
     &           inodeDatalo0:inodeDatahi0,
     &           inodeDatalo1:inodeDatahi1)
      integer icellDatalo0,icellDatalo1
      integer icellDatahi0,icellDatahi1
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0,
     &           icellDatalo1:icellDatahi1)
      integer inodeBoxlo0,inodeBoxlo1
      integer inodeBoxhi0,inodeBoxhi1
      integer i,j
      REAL*8 w,t
      w=(0.500d0)
      w=w*(0.500d0)
      do j = inodeBoxlo1,inodeBoxhi1
      do i = inodeBoxlo0,inodeBoxhi0
      t =     cellData(i-1,j-1)
     &     +  cellData(i  ,j-1)
      t = t + cellData(i-1,j  )
     &     +  cellData(i  ,j  ) 
      nodeData(i,j) = w * t
      enddo
      enddo
      return 
      end
      subroutine CELLTONODECISMVELNOSHEAR(
     &           nodeData
     &           ,inodeDatalo0,inodeDatalo1
     &           ,inodeDatahi0,inodeDatahi1
     &           ,nnodeDatacomp
     &           ,cellData
     &           ,icellDatalo0,icellDatalo1
     &           ,icellDatahi0,icellDatahi1
     &           ,nGhost
     &           ,offset
     &           ,boxlovect
     &           ,inodeBoxlo0,inodeBoxlo1
     &           ,inodeBoxhi0,inodeBoxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nnodeDatacomp
      integer inodeDatalo0,inodeDatalo1
      integer inodeDatahi0,inodeDatahi1
      REAL*8 nodeData(
     &           inodeDatalo0:inodeDatahi0,
     &           inodeDatalo1:inodeDatahi1,
     &           0:nnodeDatacomp-1)
      integer icellDatalo0,icellDatalo1
      integer icellDatahi0,icellDatahi1
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0,
     &           icellDatalo1:icellDatahi1)
      integer nGhost
      integer offset
      integer boxlovect(0:1)
      integer inodeBoxlo0,inodeBoxlo1
      integer inodeBoxhi0,inodeBoxhi1
      integer i,j
      integer nlayers, n
      integer verticalOffset
      REAL*8 w,t, nodeVal
      w=(0.500d0)
      w=w*(0.500d0)
      nlayers = (inodeDatahi0 - inodeDatalo0 +1)
      verticalOffset = -boxlovect(0)
      verticalOffset = -boxlovect(1)
      do j = inodeBoxlo1,inodeBoxhi1
      do i = inodeBoxlo0,inodeBoxhi0
      t =     cellData(i-1,j-1)
     &     +  cellData(i  ,j-1)
      t = t + cellData(i-1,j  )
     &     +  cellData(i  ,j  ) 
      nodeVal = w*t
      do n=0, nlayers-1
         nodeData(n,i,j+verticalOffset) = nodeVal
      enddo 
      enddo
      enddo
      return 
      end
