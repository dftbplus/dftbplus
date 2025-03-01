program testapp
  use fortuno_serial, only : test_list, execute_serial_cmd_app
  use test_io_indexselection, only : indexselection_tests => tests
  use test_io_tokenreader, only : tokenreader_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      indexselection_tests(),&
      tokenreader_tests()&
    ])&
  )

end program testapp
