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

module m_twostream_maxrand

      use m_data_parameters, only: ireals,iintegers,zero,one,pi
      use m_eddington, only: eddington_coeff_zdun
      use m_helper_functions, only : delta_scale_optprop
      implicit none

      private
      public delta_eddington_twostream_maxrand

    contains

      subroutine delta_eddington_twostream_maxrand(dtauf_in, w0f_in, gf_in, &
                                                   dtauc_in, w0c_in, gc_in, &
                                                   cfrac,                   &
                                                   mu0, incSolar, albedo,   &
                                                   Sf, Ednf, Eupf,          &
                                                   Sc, Ednc, Eupc) !,          &
        !                                           planckf, planckc)

        real(ireals),intent(in),dimension(:) :: dtauf_in, w0f_in, gf_in
        real(ireals),intent(in),dimension(:) :: dtauc_in, w0c_in, gc_in, cfrac

        real(ireals),intent(in) :: albedo, mu0, incSolar

        real(ireals),dimension(:),intent(out):: Sf, Ednf, Eupf 
        real(ireals),dimension(:),intent(out):: Sc, Ednc, Eupc 
        !real(ireals),dimension(:),intent(in),optional :: planckf, planckc

        real(ireals),dimension(size(dtauf_in)) :: dtauf, w0f, gf
        real(ireals),dimension(size(dtauc_in)) :: dtauc, w0c, gc
                                         
        real(ireals),dimension(size(dtauf_in)) :: a11f, a12f, a13f, a23f, a33f, g1f, g2f
        real(ireals),dimension(size(dtauc_in)) :: a11c, a12c, a13c, a23c, a33c, g1c, g2c

        real(ireals) :: b1, b2, b3, b4 ! overlap coeffs

        integer(iintegers) :: i,j,k,ke,ke1,bi
        real(ireals),allocatable :: AB (:,:)
        real(ireals),allocatable :: B (:,:)
        integer,allocatable :: IPIV(:)
        integer :: N, KLU,  KL, KU, NRHS, LDAB, LDB, INFO

        !real(ireals) :: d0,d1,c1,c2,c3 !thermal coeffs

        Sf=zero
        Ednf=zero
        Eupf=zero

        Sc=zero
        Ednc=zero
        Eupc=zero

        ke = size(dtauf)
        ke1 = ke+1
        N = 4*(ke1)
        KL = 6
        KU = 6
        NRHS = 1

        LDAB = 2*KL+KU+1
        LDB  = N
        KLU = KL+KU+1

        dtauf = dtauf_in
        w0f   = w0f_in
        gf    = gf_in

        dtauc = dtauc_in
        w0c   = w0c_in
        gc    = gc_in

        ! Delta scale cloud optical properties
        call delta_scale_optprop(dtauc, w0c, gc)

        do k=1,ke
          call eddington_coeff_zdun (dtauf(k), w0f(k),gf(k), mu0,a11f(k),a12f(k),a13f(k),a23f(k),a33f(k), g1f(k),g2f(k) )
          call eddington_coeff_zdun (dtauc(k), w0c(k),gc(k), mu0,a11c(k),a12c(k),a13c(k),a23c(k),a33c(k), g1c(k),g2c(k) )
        enddo

        ! Compute Direct radiation before hand
        Sf(1) = incSolar * (one - cfrac(1)) 
        Sc(1) = incSolar * cfrac(1)

        do k=2,ke
          
          if(cfrac(k-1).lt.one) then
            b1 = (one - max(cfrac(k), cfrac(k-1))) / (one - cfrac(k-1))
          else
            b1 = one
          endif

          if(cfrac(k-1).gt.zero) then
            b3 = min(cfrac(k), cfrac(k-1)) / cfrac(k-1)
          else
            b3 = one
          endif

          Sf(k) = a33f(k-1) * b1 * Sf(k-1) + a33c(k-1) * (one - b3) * Sc(k-1)
          Sc(k) = a33f(k-1) * (one-b1) * Sf(k-1) + a33c(k-1) * b3 * Sc(k-1)
        enddo

        ! Last level just transmit radiation through the layer til the surface...
        Sf(ke1) = a33f(ke) * Sf(ke)
        Sc(ke1) = a33c(ke) * Sc(ke)

        ! Create Space for Lapack band matrix
        allocate( IPIV(N) )
        allocate( AB (LDAB,N) )
        allocate( B (LDB,NRHS)   )
        AB   = zero
        B    = zero
        IPIV = 0

        ! Setup solar src vector
        do k=1,ke
          B(4*k-3,1) = a13f(k) * Sf(k+1) / a33f(k) ! Eupf i-1
          B(4*k+2,1) = a23f(k) * Sf(k+1) / a33f(k) ! Ednf i

          B(4*k-1,1) = a13c(k) * Sc(k+1) / a33c(k) ! Eupc i-1
          B(4*k+4,1) = a23c(k) * Sc(k+1) / a33c(k) ! Ednc i
        enddo

        B(2,1) = zero ! no Ednf at TOA
        B(4,1) = zero ! no Ednc at TOA
        B(4*ke1-3,1) = Sf(ke1) * albedo 
        B(4*ke1-1,1) = Sc(ke1) * albedo 

!       ! TODO: thermal radiation not yet implemented...
!        ! Setup thermal src vector
!        if(present(planck) ) then
!          do k=1,ke
!            if( dtau(k).gt.0.01_ireals ) then
!              d0 = planck(k)
!              d1 = ( planck(k)-planck(k+1) ) / dtau(k)
!            else
!              d0 = .5_ireals*(planck(k)+planck(k+1))
!              d1 = zero
!            endif
!            c1 = g1(k) * (d0 + d1*dtau(k))
!            c2 = g2(k) * d1
!            c3 = g1(k) * d0
!            B(2*k-1,1) = B(2*k-1,1) + ( - a11(k)*(c1+c2) - a12(k)*(c3-c2) + c2 + c3 )*pi 
!            B(2*k+2,1) = B(2*k+2,1) + ( - a12(k)*(c1+c2) - a11(k)*(c3-c2) + c1 - c2 )*pi 
!          enddo
!          ! B(2,1) = B(2,1) + zero ! no Edn at TOA
!          B(2*ke1-1,1) = B(2*ke1-1,1) + planck(ke1)*(one-albedo)*pi 
!        endif

        !diagonal entries
        do i=1,N
          j=i
          bi= KLU+i-j 
          AB( bi,j ) = one
        enddo

        ! Upward source terms....
        do k=1,ke-1
          if(cfrac(k+1).lt.one) then
            b2 = (one - max(cfrac(k), cfrac(k+1))) / (one - cfrac(k+1))
          else
            b2 = one
          endif
          if(cfrac(k+1).gt.zero) then
            b4 = min(cfrac(k), cfrac(k+1)) / cfrac(k+1)
          else
            b4 = one
          endif

        ! setting Transmission for Eupf coeffs
          i=4*k-3 ; j=i+4 ; bi= KLU+i-j ; AB(bi,j) = - (a11f(k) * b2)
          i=4*k-3 ; j=i+6 ; bi= KLU+i-j ; AB(bi,j) = - (a11f(k) * (one-b4))

        ! setting Transmission for Eupc coeffs
          i=4*k-1 ; j=i+2 ; bi= KLU+i-j ; AB(bi,j) = - (a11c(k) * (one-b2))
          i=4*k-1 ; j=i+4 ; bi= KLU+i-j ; AB(bi,j) = - (a11c(k) * b4)

        ! setting Reflection for Ednf coeffs
          i=4*k+2 ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = - (a12f(k) * b2)
          i=4*k+2 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = - (a12f(k) * (one-b4))

        ! setting Reflection for Ednc coeffs
          i=4*k+4 ; j=i-3 ; bi= KLU+i-j ; AB(bi,j) = - (a12c(k) * (one-b2))
          i=4*k+4 ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = - (a12c(k) * b4)
        enddo

        ! Downward source terms....
        do k=2,ke
          if(cfrac(k-1).lt.one) then
            b1 = (one - max(cfrac(k), cfrac(k-1))) / (one - cfrac(k-1))
          else
            b1 = one
          endif

          if(cfrac(k-1).gt.zero) then
            b3 = min(cfrac(k), cfrac(k-1)) / cfrac(k-1)
          else
            b3 = one
          endif

        ! setting Reflection for Eupf coeffs
          i=4*k-3 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = - (a12f(k) * b1)
          i=4*k-3 ; j=i+3 ; bi= KLU+i-j ; AB(bi,j) = - (a12f(k) * (one-b3))

        ! setting Reflection for Eupc coeffs
          i=4*k-1 ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = - (a12c(k) * (one-b1))
          i=4*k-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = - (a12c(k) * b3)

        ! setting Transmission for Ednf coeffs
          i=4*k+2 ; j=i-4 ; bi= KLU+i-j ; AB(bi,j) = - (a11f(k) * b1)
          i=4*k+2 ; j=i-2 ; bi= KLU+i-j ; AB(bi,j) = - (a11f(k) * (one-b3))

        ! setting Transmission for Ednc coeffs
          i=4*k+4 ; j=i-6 ; bi= KLU+i-j ; AB(bi,j) = - (a11c(k) * (one-b1))
          i=4*k+4 ; j=i-4 ; bi= KLU+i-j ; AB(bi,j) = - (a11c(k) * b3)
        enddo


        ! Last layer has no overlap rules:
          ! Transmissions...
          i=4*ke-3 ; j=i+4 ; bi= KLU+i-j ; AB(bi,j) = - a11f(ke) ! Eupf
          i=4*ke-1 ; j=i+4 ; bi= KLU+i-j ; AB(bi,j) = - a11c(ke) ! Eupc

          i=4*ke+2 ; j=i-4 ; bi= KLU+i-j ; AB(bi,j) = - a11f(ke) ! Ednf
          i=4*ke+4 ; j=i-4 ; bi= KLU+i-j ; AB(bi,j) = - a11c(ke) ! Ednc
          ! Reflections...  
          i=4*ke-3 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = - a12f(ke) ! Eupf
          i=4*ke-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = - a12c(ke) ! Eupc

          i=4*ke+2 ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = - a12f(ke) ! Ednf
          i=4*ke+4 ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = - a12c(ke) ! Ednc


        ! Surface albedo
        i=4*ke1-3 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = -albedo ! Eupf at surface is Ednf*albedo
        i=4*ke1-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = -albedo ! Eupc at surface is Ednc*albedo

        INFO=-1
        if(kind(one).eq.kind(real(one)) ) then !single_precision
          call SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
        else if(kind(one).eq.kind(dble(one)) ) then !double_precision
          call DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
        else
          print *,'Dont know which LAPACK routine to call for real kind',kind(one)
          call exit(-5)
        endif

        if(INFO.ne.0) then
          print *,'INFO',INFO
          stop 'Error in twostream calculation - lapack returned Error'
        endif

        ! retrieve result from solver
        do k=1,ke1
          Eupf(k) = B(4*k-3,NRHS) 
          Ednf(k) = B(4*k-2,NRHS) 
          Eupc(k) = B(4*k-1,NRHS) 
          Ednc(k) = B(4*k  ,NRHS) 
          if(any(isnan( [Eupf(k), Ednf(k), Eupc(k), Ednc(k)] )) ) &
              print *,'setting value for Eup,Edn',k,' LAPACK entries',B(4*k-3:4*k,1),'Eup/dn', Eupf(k), Ednf(k), Eupc(k), Ednc(k),'IPIV',IPIV(2*k-1:2*k)
        enddo

        end subroutine


end module
