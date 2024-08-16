program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_math_matrixops, only : matrixops_test_items
  implicit none

  call execute_serial_cmd_app(&
    testitems=[&
      matrixops_test_items()&
    ]&
  )

end program testapp
