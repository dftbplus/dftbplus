program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_common_atomicmass, only : atomicmass_tests => tests
  use test_common_file, only : file_tests => tests
  use test_common_memman, only : memman_tests => tests
  use test_common_schedule, only : schedule_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      atomicmass_tests(),&
      file_tests(),&
      memman_tests(),&
      schedule_tests()&
    ])&
  )

end program testapp
