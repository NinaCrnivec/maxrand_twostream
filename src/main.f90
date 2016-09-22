program main

  use m_data_parameters, only : iintegers, ireals, zero, one
  use m_twostream, only : delta_eddington_twostream
  use m_twostream_maxrand, only : delta_eddington_twostream_maxrand
  use m_nina_maxrand_wrapper, only : nina_maxrand_wrapper


  implicit none


  integer(iintegers),parameter :: ke = 60
  integer(iintegers),parameter :: ke1 = ke+1

  real(ireals),parameter :: albedo=.1, mu0=.5, incSolar=1000

  real(ireals),dimension(ke) :: cfrac
  real(ireals),dimension(ke) :: tauf, w0f, gf
  real(ireals),dimension(ke) :: tauc, w0c, gc
  real(ireals),dimension(ke) :: taupp, w0pp, gpp
  real(ireals),dimension(ke) :: kabsf, kscaf
  real(ireals),dimension(ke) :: kabsc, kscac

  real(ireals),dimension(ke1) :: Sc, Ednc, Eupc
  real(ireals),dimension(ke1) :: Sf, Ednf, Eupf

  integer(iintegers) :: k

  real :: start,finish, best_a, best_b
  integer(iintegers),parameter :: iter=1

  ! Clearsky optical properties
  tauf = one/ke*mu0
  w0f = .5_ireals
  gf  = .5_ireals
  kabsf = tauf * (one - w0f)
  kscaf = tauf * w0f

  cfrac = .3_ireals

  ! Cloudy optical properties
  tauc = tauf
  w0c = w0f
  gc  = gf
  kabsc = tauc * (one - w0c)
  kscac = tauc * w0c

  ! Plane parallel mixture of optical props
  taupp = tauf * (one - cfrac) + tauc*cfrac
  w0pp = ( kscaf + kscac ) / (tauf + tauc)
  gpp  = (kscaf*gf + kscac*gc) / (tauf + tauc) / w0pp

  ! Compute regular twostream fluxes with plane parallel approx.
  call cpu_time(start)
  do k=1,iter
    call delta_eddington_twostream(taupp, w0pp, gpp, mu0, incSolar, albedo, Sf, Ednf, Eupf)
  enddo
  call cpu_time(finish)

  print *,'twostr benchmark :: ',finish-start

  do k=1,ke1
    print *,k,'::',Sf(k), Ednf(k), Eupf(k)
  enddo

  !print *,'MaxRand computations...'
  !call cpu_time(start)
  !do k=1,iter
  !  call delta_eddington_twostream_maxrand( &
  !                    tauf, w0f, gf,        &
  !                    tauc, w0c, gc,        &
  !                    cfrac,                &
  !                    mu0, incSolar, albedo,&
  !                    Sf, Ednf, Eupf,       &
  !                    Sc, Ednc, Eupc)
  !enddo
  !call cpu_time(finish)

  !print *,'maxrand benchmark :: ',finish-start

  !do k=1,ke1
  !  print *,k,'::',Sf(k), Ednf(k), Eupf(k), '::', Sc(k), Ednc(k), Eupc(k)
  !enddo

  print *,'Nina`s MaxRand computations...'
  call cpu_time(start)
  do k=1,iter
    call nina_maxrand_wrapper( &
                      tauf, w0f, gf,        &
                      tauc, w0c, gc,        &
                      cfrac,                &
                      mu0, incSolar, albedo,&
                      Sf, Ednf, Eupf,       &
                      Sc, Ednc, Eupc)
  enddo
  call cpu_time(finish)

  print *,'nina_maxrand benchmark :: ',finish-start

  do k=1,ke1
    print *,k,'::',Sf(k), Ednf(k), Eupf(k), '::', Sc(k), Ednc(k), Eupc(k)
  enddo

end program
