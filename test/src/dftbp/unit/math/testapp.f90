program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_math_matrixops, only : matrixops_tests => tests
  implicit none

  call execute_serial_cmd_app(matrixops_tests())

end program testapp
