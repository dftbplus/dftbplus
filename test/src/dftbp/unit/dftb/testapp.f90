program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_dftb_periodic, only : periodic_test_items
  implicit none

  call execute_serial_cmd_app(&
    testitems=[&
      periodic_test_items()&
    ]&
  )

end program testapp
