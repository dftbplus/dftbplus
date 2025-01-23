program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_include_allocatablelist, only : allocatablelist_tests => tests
  use test_include_pointerlist, only : pointerlist_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      allocatablelist_tests(),&
      pointerlist_tests()&
    ])&
  )

end program testapp
