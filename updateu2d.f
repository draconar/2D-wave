c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/patchdata/fortran/pdat_m4arrdim2d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1704 $
c  Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c

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
      double precision dx(0:2-1), dt
      double precision 
     &  ucur(ifirst0-ngc0:ilast0+1+ngc0,
     &          ifirst1-ngc1:ilast1+1+ngc1),
     &  uold(ifirst0:ilast0+1,
     &          ifirst1:ilast1+1),
     &  f(ifirst0:ilast0+1,
     &          ifirst1:ilast1+1)
c output arrays:
      double precision 
     &  unew(ifirst0:ilast0+1,
     &          ifirst1:ilast1+1)
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
