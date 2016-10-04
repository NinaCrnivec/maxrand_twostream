program main

  use m_data_parameters, only : iintegers, ireals, zero, one
  use m_twostream, only : delta_eddington_twostream
  use m_twostream_maxrand, only : delta_eddington_twostream_maxrand
  use m_nina_maxrand_wrapper, only : nina_maxrand_wrapper


  implicit none


  !integer(iintegers),parameter :: ke = 4
  integer(iintegers),parameter :: ke = 5
  integer(iintegers),parameter :: ke1 = ke+1

  real(ireals),parameter :: albedo=0.1, mu0=1, incSolar=1000.0

  real(ireals),dimension(ke) :: cfrac
  real(ireals),dimension(ke) :: tauf, w0f, gf
  real(ireals),dimension(ke) :: tauc, w0c, gc
  real(ireals),dimension(ke) :: taupp, w0pp, gpp
  real(ireals),dimension(ke) :: kabsf, kscaf
  real(ireals),dimension(ke) :: kabsc, kscac

  real(ireals),dimension(ke1) :: Sc, Ednc, Eupc
  real(ireals),dimension(ke1) :: Sf, Ednf, Eupf

  integer(iintegers) :: k

  real :: start,finish
  !real :: best_a, best_b
  integer(iintegers),parameter :: iter=1

  ! Clearsky optical properties
  tauf = one/ke*mu0
  !call random_number(tauf)
  !tauf = 1e-8
  !tauf = 1e-8/ke*mu0
  !tauf(2) = 1
  w0f = .5_ireals
  gf  = .5_ireals
  kabsf = tauf * (one - w0f)
  kscaf = tauf * w0f

  cfrac = zero
  cfrac(3) = .5_ireals
  !cfrac(1) = .5_ireals
  !Or simply:
  !cfrac(1) = 0.5
  !call random_number(cfrac)
  
  
  ! TESTING p1, p2, p3, p4:
  ! 5 layers:
  !cfrac(1) = .4_ireals
  !cfrac(2) = .6_ireals
  !cfrac(3) = .3_ireals
  !cfrac(4) = .0_ireals
  !cfrac(5) = .8_ireals
  

  ! Cloudy optical properties
  tauc = zero
  w0c = zero
  gc = zero
  
  !tauc(3) = one*mu0
  !w0c(3) = .01
  !gc(3) = .9
  
  !tauc = 1
  tauc = tauf
  w0c = w0f
  gc  = gf
  kabsc = tauc * (one - w0c)
  kscac = tauc * w0c
  

  ! Plane parallel mixture of optical props
  taupp = tauf * (one - cfrac) + tauc*cfrac
  w0pp = ( kscaf + kscac ) / (tauf + tauc)
  gpp  = (kscaf*gf + kscac*gc) / (tauf + tauc) / w0pp

  
!#############
!#############
  
  !print *, ''
  !print *,'MaxRand computations...'
  
  !print *, ''
  !print *,'Nina`s MaxRand computations...'
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

  print *, ''
  print *,'nina_maxrand benchmark :: ',finish-start

  81 format(i2, 6(f12.2))
  write(6,*) "k, Edirc(k), Edirf(k), Ednc(k), Ednf(k), Eupc(k), Eupf(k)"
  do k=1,ke1
    !print *,k,'::',Sf(k), Ednf(k), Eupf(k), '::', Sc(k), Ednc(k), Eupc(k)
    write (6,81) k, Sc(k), Sf(k), Ednc(k), Ednf(k), Eupc(k), Eupf(k) 
  enddo
  
  
  
  82 format(i2, 3(f12.2))
  
  print *, ''
  write(6,*) "k, Edir(k), Edn(k), Eup(k)"
  do k=1,ke1
    !print *,k,'::',Sf(k)+Sc(k), Ednf(k)+ Ednc(k), Eupf(k)+Eupc(k)
    write(6,82) k, Sf(k)+Sc(k), Ednf(k)+ Ednc(k), Eupf(k)+Eupc(k)
  enddo
  
  
  
!#############
!#############
  
  
  
  call cpu_time(start)
  do k=1,iter
    call delta_eddington_twostream_maxrand( &
                      tauf, w0f, gf,        &
                      tauc, w0c, gc,        &
                      cfrac,                &
                      mu0, incSolar*mu0, albedo,&
                      Sf, Ednf, Eupf,       &
                      Sc, Ednc, Eupc)
  enddo
  call cpu_time(finish)

  print *, ''
  print *,''
  print *,'maxrand benchmark :: ',finish-start

  
  
  write(6,*) "k, Edirc(k), Edirf(k), Ednc(k), Ednf(k), Eupc(k), Eupf(k)"
  do k=1,ke1
    !print *,k,'::',Sf(k), Ednf(k), Eupf(k), '::', Sc(k), Ednc(k), Eupc(k)
    write (6,81) k, Sc(k), Sf(k), Ednc(k), Ednf(k), Eupc(k), Eupf(k) 
  enddo

  print *, ''
  write(6,*) "k, Edir(k), Edn(k), Eup(k)"
  do k=1,ke1
    !print *,k,'::',Sf(k)+Sc(k), Ednf(k)+ Ednc(k), Eupf(k)+Eupc(k)
    write(6,82) k, Sf(k)+Sc(k), Ednf(k)+ Ednc(k), Eupf(k)+Eupc(k)
  enddo
  
  
  
  
!#############
!#############
  
  
! Compute regular twostream fluxes with plane parallel approx.
  call cpu_time(start)
  do k=1,iter
    call delta_eddington_twostream(taupp, w0pp, gpp, mu0, incSolar*mu0, albedo, Sf, Ednf, Eupf)
    !call delta_eddington_twostream(tauf, w0f, gf, mu0, incSolar*mu0, albedo, Sf, Ednf, Eupf)
  enddo
  call cpu_time(finish)

  print *, ""
  print *,'twostr benchmark :: ',finish-start

  write(6,*) "k, Edir(k), Edn(k), Eup(k)"
  do k=1,ke1
    !print *,k,'::',Sf(k), Ednf(k), Eupf(k)
    write(6,82) k, Sf(k), Ednf(k), Eupf(k)
  enddo
  
  
!#############
!#############
  
end program
