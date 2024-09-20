program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_include_allocatablelist, only : allocatablelist_test_items
  use test_include_pointerlist, only : pointerlist_test_items
  implicit none

  call execute_serial_cmd_app(&
    testitems=[&
      allocatablelist_test_items(),&
      pointerlist_test_items()&
    ]&
  )

end program testapp
