program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_wavegrid_basis_lut, only : lut_tests => tests
  use test_wavegrid_basis_slater, only : sto_tests => tests
  use test_wavegrid_basis_spharmonics, only : ylm_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      lut_tests(), sto_tests(), ylm_tests()&
  ])&
  )

end program testapp
