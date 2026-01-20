program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_math_matrixops, only : matrixops_tests => tests
  use test_math_interpolation, only : interpolation_tests => tests
  implicit none

  call execute_serial_cmd_app(interpolation_tests())
  call execute_serial_cmd_app(matrixops_tests())

end program testapp
