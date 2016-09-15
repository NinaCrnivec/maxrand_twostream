module m_data_parameters
  double precision, parameter :: dummy_double=0
  integer,parameter :: dummy_integer=0

  integer,parameter :: ireals=kind(dummy_double),    &
                       iintegers=kind(dummy_integer), &
                       ireal128 = selected_real_kind(33, 4931)

  real(ireals),parameter :: zero=0, one=1, pi=3.141592653589793_ireals
end module

