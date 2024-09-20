program testapp
  use fortuno_serial, only : execute_serial_cmd_app
  use test_common_atomicmass, only : atomicmass_test_items
  use test_common_file, only : file_test_items
  use test_common_memman, only : memman_test_items
  use test_common_schedule, only : schedule_test_items
  implicit none

  call execute_serial_cmd_app(&
    testitems=[&
      atomicmass_test_items(),&
      file_test_items(),&
      memman_test_items(),&
      schedule_test_items()&
    ]&
  )

end program testapp
