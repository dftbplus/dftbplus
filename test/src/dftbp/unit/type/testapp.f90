program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_type_typegeometryhsd, only : typegeometryhsd_tests => tests
  implicit none

  call execute_serial_cmd_app(typegeometryhsd_tests())

end program testapp
