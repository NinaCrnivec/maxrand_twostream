module m_nina_maxrand_wrapper

  use iso_c_binding

  use m_data_parameters, only: ireals,iintegers,zero,one,pi
  implicit none

  private
  public nina_maxrand_wrapper

  interface
    integer(c_int) function nina_maxrand( nlyr, &
        dtauf, w0f, gf, &
        dtauc, w0c, gc, &
        cfrac,                   &
        mu0, incSolar, albedo,   &
        Sf, Ednf, Eupf,          &
        Sc, Ednc, Eupc) bind ( c )

      use iso_c_binding
      integer ( c_int ), VALUE :: nlyr
      real ( c_double ) :: dtauf(*), w0f(*), gf(*)
      real ( c_double ) :: dtauc(*), w0c(*), gc(*)
      real ( c_double ) :: cfrac(*)
      real ( c_double ), VALUE :: mu0, incSolar, albedo
      real ( c_double ) :: Sf(*), Ednf(*), Eupf(*)
      real ( c_double ) :: Sc(*), Ednc(*), Eupc(*)

    end function
  end interface

contains

  subroutine nina_maxrand_wrapper( &
      dtauf, w0f, gf, &
      dtauc, w0c, gc, &
      cfrac,                   &
      mu0, incSolar, albedo,   &
      Sf, Ednf, Eupf,          &
      Sc, Ednc, Eupc)

    real(ireals),intent(in),dimension(:) :: dtauf, w0f, gf
    real(ireals),intent(in),dimension(:) :: dtauc, w0c, gc, cfrac

    real(ireals),intent(in) :: albedo, mu0, incSolar

    real(ireals),dimension(:),intent(out):: Sf, Ednf, Eupf 
    real(ireals),dimension(:),intent(out):: Sc, Ednc, Eupc 

    integer(c_int) :: nlyr, ierr

    nlyr = size(dtauf)

    ierr = nina_maxrand( nlyr, &
      dtauf, w0f, gf, &
      dtauc, w0c, gc, &
      cfrac,                   &
      mu0, incSolar, albedo,   &
      Sf, Ednf, Eupf,          &
      Sc, Ednc, Eupc)


  end subroutine


end module
