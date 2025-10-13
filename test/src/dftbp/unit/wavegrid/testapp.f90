program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_wavegrid_simple, only : simple_tests => tests
  use test_wavegrid_radial, only : radial_tests => tests
  use test_wavegrid_spharmonics, only : spharm_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      simple_tests(), radial_tests(), spharm_tests() &
    ])&
  )

end program testapp
