define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine updateu(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngc0,ngc1,
     &  dx, dt,
     &  unew,
     &  ucur,
     &  uold,
     &  f)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ngc0,ngc1
      double precision dx(0:NDIM-1), dt
      double precision 
     &  ucur(NODE2dVECG(ifirst,ilast,ngc)),
     &  uold(NODE2d(ifirst,ilast,0)),
     &  f(NODE2d(ifirst,ilast,0))
c output arrays:
      double precision 
     &  unew(NODE2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i,j
      double precision xdiff, ydiff

c
c  Computes updated u for 2nd order acc linear wave eqn:
c
c     unew(i) = dt^2/dx^2*(ucur(i+1) - 2*ucur(i) + ucur(i-1)) 
c                 + 2*ucur(i) - uold(i) + dt^2*f
c
c  Note that ucur requires ghosts while unew and uold do not.  
c 
      do j=ifirst1,ilast1+1
         do i=ifirst0,ilast0+1

c........compute central differences in X and Y

            xdiff = (ucur(i+1,j)-2.0*ucur(i,j)+ucur(i-1,j))/dx(0)**2
            ydiff = (ucur(i,j+1)-2.0*ucur(i,j)+ucur(i,j-1))/dx(1)**2

c........compute updated u

c  1d
c            unew(i,j) = dt**2 * xdiff + 
c     &                  2.0*ucur(i,j) - uold(i,j) + dt**2 * f(i,j)
c  2d
            unew(i,j) = dt**2 * (xdiff + ydiff) + 
     &                  2.0*ucur(i,j) - uold(i,j) + dt**2 * f(i,j)

         enddo
      enddo
c
      return
      end
c
c
