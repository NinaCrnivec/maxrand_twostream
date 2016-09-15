!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_helper_functions
      use m_data_parameters,only : iintegers,ireals,pi,zero,one

      implicit none

    contains
      pure elemental subroutine swap(x,y)
          real(ireals),intent(inout) :: x,y
          real(ireals) :: tmp
          tmp = x
          x = y
          y = tmp
      end subroutine
      pure elemental subroutine inc(x,i)
          real(ireals),intent(inout) :: x
          real(ireals),intent(in) :: i
          x=x+i
      end subroutine

      pure function gradient(v)
          real(ireals),intent(in) :: v(:)
          real(ireals) :: gradient(size(v)-1)
          gradient = v(2:size(v))-v(1:size(v)-1)
      end function

      pure function meanvec(v)
          real(ireals),intent(in) :: v(:)
          real(ireals) :: meanvec(size(v)-1)
          meanvec = (v(2:size(v))+v(1:size(v)-1))*.5_ireals
      end function

      pure function norm(v)
          real(ireals) :: norm
          real(ireals),intent(in) :: v(:)
          norm = sqrt(dot_product(v,v))
      end function

      elemental function deg2rad(deg)
          real(ireals) :: deg2rad
          real(ireals),intent(in) :: deg
          deg2rad = deg *pi/180._ireals
      end function
      elemental function rad2deg(rad)
          real(ireals) :: rad2deg
          real(ireals),intent(in) :: rad
          rad2deg = rad /pi*180._ireals
      end function

      pure function rmse(a,b)
          real(ireals) :: rmse(2)
          real(ireals),intent(in) :: a(:),b(:)
          rmse(1) = sqrt( mean( (a-b)**2 ) )
          rmse(2) = rmse(1)/max( mean(b), epsilon(rmse) )
      end function

      pure function mean(arr)
          real(ireals) :: mean
          real(ireals),intent(in) :: arr(:)
          mean = sum(arr)/size(arr)
      end function

      elemental logical function approx(a,b,precision)
          real(ireals),intent(in) :: a,b
          real(ireals),intent(in),optional :: precision
          real(ireals) :: factor
          if(present(precision) ) then
            factor = precision
          else
            factor = 10._ireals*epsilon(b)
          endif
          if( a.le.b+factor .and. a.ge.b-factor ) then
            approx = .True.
          else
            approx = .False.
          endif
      end function
      elemental logical function rel_approx(a,b,precision)
          real(ireals),intent(in) :: a,b
          real(ireals),intent(in),optional :: precision
          real(ireals) :: factor,rel_error
          if(present(precision) ) then
            factor = precision
          else
            factor = 10*epsilon(b)
          endif
          rel_error = abs( (a-b)/ max(epsilon(a), ( (a+b)*.5_ireals ) ) )

          if( rel_error .lt. precision ) then
            rel_approx = .True.
          else
            rel_approx = .False.
          endif
      end function

      elemental subroutine delta_scale( kabs,ksca,g,factor ) 
          real(ireals),intent(inout) :: kabs,ksca,g ! kabs, ksca, g
          real(ireals),intent(in),optional :: factor
          real(ireals) :: dtau, w0
          dtau = max( kabs+ksca, epsilon(dtau) )
          w0   = ksca/dtau
          g    = g

          if(present(factor)) then
            call delta_scale_optprop( dtau, w0, g, factor)
          else
            call delta_scale_optprop( dtau, w0, g)
          endif

          kabs= dtau * (one-w0)
          ksca= dtau * w0
      end subroutine
      elemental subroutine delta_scale_optprop( dtau, w0, g,factor) 
          real(ireals),intent(inout) :: dtau,w0,g
          real(ireals),intent(in),optional :: factor
          real(ireals) :: f

          g = min( g, one-epsilon(g)*10)
          if(present(factor)) then
            f = factor
          else
            f = g**2
          endif
          dtau = dtau * ( one - w0 * f )
          g    = ( g - f ) / ( one - f )
          w0   = w0 * ( one - f ) / ( one - f * w0 )
      end subroutine

      function cumsum(arr)
          real(ireals),intent(in) :: arr(:)
          real(ireals) :: cumsum(size(arr))
          integer :: i
          cumsum(1) = arr(1)
          do i=2,size(arr)
            cumsum(i) = cumsum(i-1) + arr(i)
          enddo
      end function

      subroutine read_ascii_file_2d(filename, arr, ncolumns, skiplines, ierr)
          character(len=*),intent(in) :: filename
          integer(iintegers),intent(in) :: ncolumns
          integer(iintegers),intent(in),optional :: skiplines

          real(ireals),allocatable,intent(out) :: arr(:,:)

          integer(iintegers) :: ierr

          real :: line(ncolumns)

          integer(iintegers) :: unit, nlines, i, io
          logical :: file_exists=.False.

          ierr=0
          inquire(file=filename, exist=file_exists)

          if(.not.file_exists) then
            print *,'File ',trim(filename), 'does not exist!'
            ierr=1
            return
          endif

          open(newunit=unit, file=filename)
          if(present(skiplines)) then
              do i=1,skiplines
                  read(unit,*)
              enddo
          endif

          nlines = 0
          do
              read(unit, *, iostat=io) line
              !print *,'line',line
              if (io/=0) exit
              nlines = nlines + 1
          end do

          rewind(unit)
          if(present(skiplines)) then
              do i=1,skiplines
                  read(unit,*)
              enddo
          endif

          allocate(arr(nlines,ncolumns))

          do i=1,nlines
              read(unit, *, iostat=io) line
              arr(i,:) = line
          end do

          close(unit)
          print *,'I read ',nlines,'lines'
      end subroutine

      function search_sorted_bisection(arr,val) ! return index+residula i where arr(i) .gt. val
        real(ireals) :: search_sorted_bisection
        real(ireals),intent(in) :: arr(:)
        real(ireals),intent(in) :: val
        real(ireals) :: loc_increment
        integer(iintegers) :: i,j,k

        i=lbound(arr,1)
        j=ubound(arr,1)

        do
          k=(i+j)/2
          if (val < arr(k)) then
            j=k
          else
            i=k
          endif
          if (i+1 >= j) then ! only single or tuple left
            ! i is left bound and j is right bound index
            if(i.eq.j) then
              loc_increment = zero
            else
              loc_increment = (val - arr(i)) / ( arr(j) - arr(i) )
            endif
            search_sorted_bisection= min(max(one*lbound(arr,1), i + loc_increment), one*ubound(arr,1)) ! return `real-numbered` location of val
            exit
          endif
        end do
      end function

  end module
