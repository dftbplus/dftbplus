program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_dftb_periodic, only : periodic_tests => tests
  implicit none

  call execute_serial_cmd_app(periodic_tests())

end program testapp
