program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_math_matrixops, only : matrixops_tests => tests
  use test_math_lapackroutines, only : lapack_svd_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      matrixops_tests(),&
      lapack_svd_tests()&
    ])&
  )

end program testapp
