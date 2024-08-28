#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine AVERAGEHO(
     &           coarse
     &           ,icoarselo0,icoarselo1
     &           ,icoarsehi0,icoarsehi1
     &           ,fine
     &           ,ifinelo0,ifinelo1
     &           ,ifinehi0,ifinehi1
     &           ,lap
     &           ,ilaplo0,ilaplo1
     &           ,ilaphi0,ilaphi1
     &           ,iblo0,iblo1
     &           ,ibhi0,ibhi1
     &           ,ilapBoxlo0,ilapBoxlo1
     &           ,ilapBoxhi0,ilapBoxhi1
     &           ,nref
     &           ,ifineRefBoxlo0,ifineRefBoxlo1
     &           ,ifineRefBoxhi0,ifineRefBoxhi1
     &           ,doHO
     &           ,doAverage
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1)
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1)
      integer ilaplo0,ilaplo1
      integer ilaphi0,ilaphi1
      REAL_T lap(
     &           ilaplo0:ilaphi0,
     &           ilaplo1:ilaphi1)
      integer iblo0,iblo1
      integer ibhi0,ibhi1
      integer ilapBoxlo0,ilapBoxlo1
      integer ilapBoxhi0,ilapBoxhi1
      integer nref
      integer ifineRefBoxlo0,ifineRefBoxlo1
      integer ifineRefBoxhi0,ifineRefBoxhi1
      integer doHO
      integer doAverage
      integer var
      integer ic0,ic1
      integer ifine0,ifine1
      integer ifinecell0,ifinecell1
      integer dir
      REAL_T weight
#if (CH_SPACEDIM > 1)
#endif
#if (CH_SPACEDIM > 2)
#endif
      
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

         coarse(ic0,ic1) = zero
      
      enddo
      enddo
      
      do ifine1 = ilapBoxlo1,ilapBoxhi1
      do ifine0 = ilapBoxlo0,ilapBoxhi0

        if ((CH_SPACEDIM.eq.2).and.(doHO.ne.0)) then
#if (CH_SPACEDIM == 2)
           weight = (nref+1)*nref/eight/four
           lap(ifine0,ifine1) = -weight
     &                               *(-four*fine(ifine0,ifine1)
     &                                +fine(ifine0+1,ifine1)
     &                                +fine(ifine0-1,ifine1)
     &                                +fine(ifine0,ifine1+1)
     &                                +fine(ifine0,ifine1-1))
#endif
        else
          lap(ifine0,ifine1) = 0.0
        endif
      
      enddo
      enddo
      if (doAverage .eq.1) then
         weight = one/(nref**CH_SPACEDIM)
      else
         weight = one
      endif
      
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

        
      do ifinecell1 = ifineRefBoxlo1,ifineRefBoxhi1
      do ifinecell0 = ifineRefBoxlo0,ifineRefBoxhi0

        
          ifine0 = nref*ic0 + ifinecell0
          ifine1 = nref*ic1 + ifinecell1
          coarse(ic0,ic1) = coarse(ic0,ic1)
     &                 + weight*fine(ifine0,ifine1)
     &                 + weight*lap(ifine0,ifine1)
        
      enddo
      enddo
      
      enddo
      enddo
      return
      end
