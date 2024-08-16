program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_type_typegeometryhsd, only : typegeometryhsd_test_items
  implicit none

  call execute_serial_cmd_app(&
    testitems=[&
      typegeometryhsd_test_items()&
    ]&
  )

end program testapp
