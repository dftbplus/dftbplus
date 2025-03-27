program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_mixer_anderson, only : anderson_tests => tests
  use test_mixer_broyden, only : broyden_tests => tests
  use test_mixer_diis, only : diis_tests => tests
  use test_mixer_simple, only : simple_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      anderson_tests(),&
      broyden_tests(),&
      diis_tests(),&
      simple_tests()&
    ])&
  )

end program testapp
