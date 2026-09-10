program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_elecsolvers_partialdiag, only : partialdiag_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      partialdiag_tests()&
    ])&
  )

end program testapp
