!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program makecube
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: i,j,k, nx,ny,nz, narg, ln, err
  real(dp), dimension(:,:,:), allocatable :: phi3d, phi3d_0
  real(dp), dimension(:), allocatable :: x,y,z
  real(dp), dimension(3) :: or
  real(dp), parameter :: au = 0.529177_dp
  character(256) :: filebox,filex,filey,filez,filename,refname
  logical :: refpot
  integer :: fp

  narg=command_argument_count()

  if (.not.(narg.eq.1 .or. narg.eq.3 .or. narg.eq.8)) then
    write(*,*) 'usage:'
    write(*,*) 'makecube pot_file [-r refpot] [-b boxfile xfile yfile zfile] '
    stop
  endif

  call get_command_argument(1,filename,ln,err)

  call get_command(filebox, ln, err)
  if (ln<0 .or. err>0) then
    stop 'Internal error: command line too long'
  endif

  k = index(filebox,"-r")
  if (k > 0) then
    refpot = .true.
    call get_command_argument(3,refname,ln,err)
  else
    refpot = .false.
  endif

  k = index(filebox,"-b")
  if (k.eq.0) then
    filebox="box3d.dat"
    filex="Xvector.dat"
    filey="Yvector.dat"
    filez="Zvector.dat"
  else
    call get_command_argument(narg-3,filebox,ln,err)
    call get_command_argument(narg-2,filex,ln,err)
    call get_command_argument(narg-1,filey,ln,err)
    call get_command_argument(narg,filez,ln,err)
  endif


  open(newunit=fp,file=trim(filebox))
  read(fp,*) nx,ny,nz
  close(fp)

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))
  allocate(phi3d(nx,ny,nz))
  if (refpot) allocate(phi3d_0(nx,ny,nz))

  open(newunit=fp,file=filex)
  read(fp,*) x
  close(fp)

  open(newunit=fp,file=filey)
  read(fp,*) y
  close(fp)

  open(newunit=fp,file=filez)
  read(fp,*) z
  close(fp)

  open(newunit=fp,file=filename)

  do i=1,nx
    do j=1,ny
      do k=1,nz
        read(fp,*) phi3d(i,j,k)
      enddo
    enddo
  enddo
  close(fp)

  if (refpot) then
    open(newunit=fp,file=refname)

    do i=1,nx
      do j=1,ny
        do k=1,nz
          read(fp,*) phi3d_0(i,j,k)
        enddo
      enddo
    enddo
    close(fp)
  endif

  or(1) = x(1) !(x(1) - x(nx))/2.d0
  or(2) = y(1) !(y(1) - y(ny))/2.d0
  or(3) = z(1) !(z(1) - z(nz))/2.d0

  write(*,*) 'CUBE'
  write(*,*) 'x, y, z'
  write(*,'(i4,3f12.5)') 1,or(1)/au,or(2)/au,or(3)/au
  write(*,'(i4,3f12.5)') nx,(x(2)-x(1))/au,0.0,0.0
  write(*,'(i4,3f12.5)') ny,0.0,(y(2)-y(1))/au,0.0
  write(*,'(i4,3f12.5)') nz,0.0,0.0,(z(2)-z(1))/au
  write(*,'(i1,4f12.5)') 1,0.0,0.0,0.0,0.0


  do i=1,nx
    do j=1,ny
      do k=1,nz
        if (refpot) then
          write(*,'(E13.5)',advance='NO') phi3d(i,j,k)-phi3d_0(i,j,k)
        else
          write(*,'(E13.5)',advance='NO') phi3d(i,j,k)
        endif
        if (mod(k-1,6) .eq. 5) write(*,*)
      enddo
      write(*,*)
    enddo
  enddo

  open(newunit=fp,file="poissonbox.jmol")

  write(fp,'(a,3(F10.4),a,3(F10.4),a,3(F10.4),a,3(F10.4),a)') &
      "draw p1 PLANE {",x(1),y(1),z(1),"}{", x(nx),y(1),z(1),"}{", &
      x(nx),y(ny),z(1),"}{",x(1),y(ny),z(1),"} color yellow translucent"

  write(fp,'(a,3(F10.4),a,3(F10.4),a,3(F10.4),a,3(F10.4),a)') &
      "draw p2 PLANE {",x(1),y(1),z(nz),"}{",x(nx),y(1),z(nz),"}{", &
      x(nx),y(ny),z(nz),"}{",x(1),y(ny),z(nz),"} color yellow translucent"

  write(fp,'(a,3(F10.4),a,3(F10.4),a)') "draw l1 {",x(1),y(1),z(1),"}{",x(1),y(1),z(nz),"}"
  write(fp,'(a,3(F10.4),a,3(F10.4),a)') "draw l2 {",x(nx),y(1),z(1),"}{",x(nx),y(1),z(nz),"}"
  write(fp,'(a,3(F10.4),a,3(F10.4),a)') "draw l3 {",x(nx),y(ny),z(1),"}{",x(nx),y(ny),z(nz),"}"
  write(fp,'(a,3(F10.4),a,3(F10.4),a)') "draw l4 {",x(1),y(ny),z(1),"}{",x(1),y(ny),z(nz),"}"

  close(fp)


  deallocate(x,y,z)
  deallocate(phi3d)
  if (refpot) deallocate(phi3d_0)

end program makecube
