program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_io_indexselection, only : indexselection_test_items
  use test_io_tokenreader, only : tokenreader_test_items
  implicit none

  call execute_serial_cmd_app(&
    testitems=[&
      indexselection_test_items(),&
      tokenreader_test_items()&
    ]&
  )

end program testapp
