program buildwire
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: i,d(3),pl_atm,tot_atm,num_pl,err,k
  character(100) :: gen_file
  character(2) :: period
  character(70) :: atm_spec
  character(100) :: cell_centre, tmpStr
  real(dp) :: cell(3,3)
  integer, ALLOCATABLE, DIMENSION (:) :: n_atm,typ_atm
  real(dp), ALLOCATABLE, DIMENSION (:) :: X,Y,Z
  character(100) :: arg
  logical :: do_super
  integer :: iargc

  iargc = command_argument_count()
  if (iargc < 3) then
    write(*,*) 'buildwire pl.gen dir npls [-s]'
    write(*,*) 'pl.gen: gen file with PL definition (supercell)'
    write(*,*) 'dir: 1=x, 2=y, 3=z'
    write(*,*) 'npls: Number of pls in the scattering region'
    write(*,*) '-s: (optional) makes a supercell structure'
    stop
  end if

  d = 0
  call get_command_argument(1, arg)
  read(arg,*) gen_file
  call get_command_argument(2, arg)
  read(arg,*) k
  d(k) = 1
  call get_command_argument(3, arg)
  read(arg,*) num_pl

  do_super = .false.
  if (iargc == 4) then
    call get_command_argument(4, arg)
    if (trim(arg) == "-s") then
      do_super = .true.
    end if
  end if

  open(30,file=trim(gen_file))

  read(30,*) pl_atm, period
  read(30,'(A)') atm_spec

  ALLOCATE(X(pl_atm),stat=err)
  IF (err /= 0) STOP 'no space for allocation (X)'

  ALLOCATE(Y(pl_atm),stat=err)
  IF (err /= 0) STOP 'no space for allocation (Y)'

  ALLOCATE(Z(pl_atm),stat=err)
  IF (err /= 0) STOP 'no space for allocation (Z)'

  ALLOCATE(n_atm(pl_atm),stat=err)
  IF (err /= 0) STOP 'no space for allocation (n_atm)'

  ALLOCATE(typ_atm(pl_atm),stat=err)
  IF (err /= 0) STOP 'no space for allocation (typ_atm)'


  do i = 1,pl_atm
     read(30,*) n_atm(i),typ_atm(i),X(i),Y(i),Z(i)
  end do

  read(30,'(A)') cell_centre
  do i = 1,3
     read(30,*) cell(i,1),cell(i,2),cell(i,3)
  end do

  close(30)


  tot_atm=pl_atm * (num_pl + 4)


  open(30,file='Ordered_'//trim(gen_file))

  if (.not.do_super) then
     period = "C"
  end if

  write(30,'(I5,A4)') tot_atm, period
  write(30,'(A)') trim(atm_spec)
  do k = 1,num_pl
     do i = 1,pl_atm
        write(30,'(I5,I5,F20.12,F20.12,F20.12)') n_atm(i),typ_atm(i), &
                               X(i)+(k-1)*cell(1,1)*d(1),&
                               Y(i)+(k-1)*cell(2,2)*d(2), &
                               Z(i)+(k-1)*cell(3,3)*d(3)
     end do
  enddo

  ! build I contact
  do k = num_pl+1,num_pl+2
     do i = 1,pl_atm
        write(30,'(I5,I5,F20.12,F20.12,F20.12)') n_atm(i),typ_atm(i), &
                               X(i)+(k-1)*cell(1,1)*d(1),&
                               Y(i)+(k-1)*cell(2,2)*d(2), &
                               Z(i)+(k-1)*cell(3,3)*d(3)
     end do
  enddo

  ! build II contact
  do k = 0,-1,-1
     do i = 1,pl_atm
        write(30,'(I5,I5,F20.12,F20.12,F20.12)') n_atm(i),typ_atm(i), &
                               X(i)+(k-1)*cell(1,1)*d(1),&
                               Y(i)+(k-1)*cell(2,2)*d(2), &
                               Z(i)+(k-1)*cell(3,3)*d(3)
     end do
  enddo

  if (do_super) then
    write(30,'(A)') cell_centre
    do i = 1,3
        if(d(1).eq.1) write(30,*) cell(i,1)*(num_pl+4)*d(1),cell(i,2),cell(i,3)
        if(d(2).eq.1) write(30,*) cell(i,1),cell(i,2)*(num_pl+4)*d(2),cell(i,3)
        if(d(3).eq.1) write(30,*) cell(i,1),cell(i,2),cell(i,3)*(num_pl+4)*d(3)
    end do
  end if

  close(30)

  write(*,*) 'structure built'
  write(*,*) 'Input for dftb+:'
  write(*,*)
  write(*,*) 'Transport{'
  write(*,*) '  Device{'
  write(*,"(1x,A,I0,' ',I0)") '    AtomRange = ',1, pl_atm*num_pl
  write(*,FMT='(1x,a)', advance='NO') '    FirstLayerAtoms ='
  do i=1,num_pl
    write(tmpStr,"(I0)") (i-1)*pl_atm+1
    write(*,'(1X,A)', advance='NO')trim(tmpStr)
  end do
  write(*,*)
  write(*,*) '  }'

  write(*,*) '  Contact{'
  write(*,*) '    Id = "source"'
  write(*,"(1x,A,I0,' ',I0)") '    AtomRange = ', pl_atm*num_pl+1, pl_atm*(num_pl+2)
  write(*,*) '  }'

  write(*,*) '  Contact{'
  write(*,*) '    Id = "drain"'
  write(*,"(1x,A,I0,' ',I0)") '    AtomRange = ', pl_atm*(num_pl+2)+1, pl_atm*(num_pl+4)
  write(*,*) '  }'

  write(*,*) '  Task= contactHamiltonian{'
  write(*,*) '    contactId = "source"'
  write(*,*) '  }'
  write(*,*) '}'
  write(*,*)

end program
