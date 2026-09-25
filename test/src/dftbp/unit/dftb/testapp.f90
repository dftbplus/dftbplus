program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_dftb_periodic, only : periodic_tests => tests
  use test_dftb_coulomb, only : coulomb_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      periodic_tests(),&
      coulomb_tests()&
    ])&
  )

end program testapp
