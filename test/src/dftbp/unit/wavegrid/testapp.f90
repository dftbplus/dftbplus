program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_wavegrid_spotcheck, only : spotcheck_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      spotcheck_tests() &
    ])&
  )

end program testapp
