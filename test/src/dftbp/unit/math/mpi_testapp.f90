program mpi_testapp
  use fortuno_mpi, only : execute_mpi_cmd_app
  use test_math_mpi_matrixops, only : tests
  implicit none

  call execute_mpi_cmd_app(tests())

end program mpi_testapp
