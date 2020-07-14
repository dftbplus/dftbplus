! PGI compilers having trouble with huge(1) in loops
!
! Known to fail:
!   PGI 17.10
!
! Workaround: use huge(1) - 1
!
program test
  implicit none

  integer :: ii

  print *, "Should see a large number of values"

  do ii = 0, huge(1)
    write(*,*)ii
  end do

end program test

