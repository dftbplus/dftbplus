module gprecision
  use, intrinsic :: iso_fortran_env, only : int64

  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,100)
  integer, parameter :: long = int64

  real(dp), parameter :: EPS=1.d-10

end module gprecision
