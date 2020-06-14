!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to make socket contact with an external code and
!! communicate data back and forward from DFTB+ to the external code.
module dftbp_ipisocket
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_message
  use dftbp_fsockets
  use dftbp_logger, only : LogWriter
  implicit none
  private

  public :: IpiSocketCommInp
  public :: IpiSocketComm, IpiSocketComm_init
  public :: IPI_PROTOCOLS


  !> Input for initialising IpiSocketComm.
  type :: IpiSocketCommInp

    !> Number of atoms for which data is exchanged
    integer :: nAtom

    !> Host name
    character(:), allocatable :: host

    !> Verbosity level of detail on communication.
    integer :: verbosity

    !> Protocol type of message headers and data to use (currently only
    !! IPI_PROTOCOL1 understood).
    integer :: protocol

    !> Port to connect to if using an internet protocol, if -1, its a file
    !! system connection
    integer :: port

  end type IpiSocketCommInp


  !> Communicator for i-Pi communication via sockets
  type :: IpiSocketComm
    private

    !> used to log messages
    type(LogWriter) :: logger

    !> level of verbosity
    integer :: verbosity

    !> socket number
    integer :: socket

    !> expected number of atoms
    integer :: nAtom

    !> Initialisation of variables
    logical :: tInit = .false.
  contains

    !> send data out
    procedure :: send

    !> receive data in
    procedure :: receive

    !> shut the socket down
    procedure :: shutdown
  end type IpiSocketComm


  !> Constructor for IpiSocketComm.
  interface IpiSocketComm
    module procedure construct
  end interface IpiSocketComm


  !> Enumerate possible protocols
  type :: IpiProtocolsEnum
    integer :: IPI_1
  end type IpiProtocolsEnum


  !> Available socket protocol types
  type(IpiProtocolsEnum), parameter :: IPI_PROTOCOLS =&
      & IpiProtocolsEnum(1)


  !> Length of strings expected for i-pi messages
  integer, parameter :: IPI_MSGLEN = 12

contains


  !> Construct IpiSocketComm instance.
  !!
  subroutine IpiSocketComm_init(this, input)

    !> Instance.
    type(IpiSocketComm), intent(out) :: this

    !> Input data.
    type(IpiSocketCommInp), intent(in) :: input

    logical :: tUnix
    character(lc) :: msg

    @:ASSERT(.not. this%tInit)
    @:ASSERT(input%natom > 0)
    @:ASSERT(input%verbosity >= 0)

    this%nAtom = input%nAtom
    this%logger = LogWriter(input%verbosity)

    if (input%protocol /= IPI_PROTOCOLS%IPI_1) then
      write(msg, '(A,1X,I0)') 'Unknown message protocol', input%protocol
    end if

    call this%logger%write('socketCreate: Opening socket for two-way&
        & communication with a server.', 1)

    ! If port > 0: internet port, otherwise local Unix socket
    tUnix = input%port < 1
    if (tUnix) then
      call this%logger%write('Establishing UNIX socket connection to'&
          & // trim(input%host), 1)
      call connect_unix_socket(this%socket, input%host)
    else
      call this%logger%write('Establishing an internet connection to', 1)
      call this%logger%write('Host: ' // trim(input%host), 1)
      write(msg, '(A,I0)') 'Port: ', input%port
      call this%logger%write(msg, 1)
      call connect_inet_socket(this%socket, input%host, input%port)
    end if

    call this%logger%write('socketCreate: ...Done', 1)
    this%tInit = .true.

  end subroutine IpiSocketComm_init


  !> Construct IpiSocketComm instance.
  !!
  function construct(input) result(this)

    !> Input data
    type(IpiSocketCommInp), intent(in) :: input

    !> Instance
    type(IpiSocketComm) :: this

    call IpiSocketComm_init(this, input)

  end function construct


  !> Receives data from an external i-Pi program via a socket.
  !!
  !! All data in atomic units, and currently assumes the number
  !! of atoms is the same as passed at construction/initialisation.
  !!
  subroutine receive(this, coord, cell, tStop)

    !> Instance.
    class(IpiSocketComm), intent(in) :: this

    !> Atomic coordinates.
    real(dp), intent(inout) :: coord(:,:)

    !> Cell lattice vectors.
    real(dp), intent(inout) :: cell(3, 3)

    !> Halt DFTB+
    logical, intent(out) :: tStop

    character(lc) :: msg
    character(len=IPI_MSGLEN) :: header, buffer
    integer :: nAtom

    ! single precision in the communication interface for some reason
    real(rdp), allocatable :: commsBuffer1(:)
    real(rdp) :: commsBuffer2(9)

    @:ASSERT(this%tInit)
    @:ASSERT(size(coord, dim=1) == 3)

    tStop = .false.

    nAtom = size(coord, dim=2)
    if (nAtom /= this%nAtom) then
      write(msg, '(1X,A,2I4)') 'Mismatch in number of atoms in socketRetrieve',&
          & nAtom, this%nAtom
      call error(msg)
    end if

    allocate(commsBuffer1(this%nAtom * 3))
    call this%logger%write('socketRetrieve: Retrieving data from socket... ', 1)

    ! wait for anything other than 'STATUS' state from the interface, returning state 'READY' in the
    ! meanwhile
    do while (.true.)
      call readbuffer(this%socket, header)
      call this%logger%write('ipisocket%receive: read from socket: ' // trim(header), 3)
      if (trim(header) /= 'STATUS') then
        exit
      end if

      buffer = 'READY'
      call writebuffer(this%socket, buffer)
      call this%logger%write('ipisocket%receive: write to socket: READY', 3)
    end do

    ! expecting positions data
    select case (trim(header))
    case ('POSDATA')
      ! expected return during run
    case ('EXIT')
      call warning("ipisocket%receive: EXIT. Halting DFTB+.")
      tStop = .true.
      return
    case default
      call error("ipisocket%receive: Unexpected message from server, received '" &
          & // trim(header) // "'")
    end select

    ! lattice vector data
    call readbuffer(this%socket, commsBuffer2)
    cell(:,:) = reshape(commsBuffer2, [3, 3])

    call this%logger%write('ipisocket%receive: read from socket: cell', 3)
    call this%logger%write(cell, 4, '(f12.6)')

    ! inverse lattice vectors
    call readbuffer(this%socket, commsBuffer2)
    call this%logger%write('ipisocket%receive: read from socket: inverse of cell', 3)

    ! number of atomic coordinates
    call readbuffer(this%socket, natom)
    call this%logger%write('ipisocket%receive: read from socket: number of atoms', 3)
    write(msg, '(I0)') nAtom
    call this%logger%write(msg, 4)

    if (nAtom /= this%nAtom) then
      write(msg, '(1X,A,2I4)')'Mismatch in number of atoms received', nAtom,this%nAtom
      call error(msg)
    end if

    ! read actual coordinates
    call readbuffer(this%socket, commsBuffer1)
    coord = reshape(commsBuffer1, [3, natom])
    call this%logger%write('ipisocket%receive: read from socket: atomic positions', 3)
    call this%logger%write(coord, 5, '(f12.6)')
    call this%logger%write('ipisocket%receive: Done', 1)

  end subroutine receive


  !> Send data to an external program via a socket
  !!
  !! All data in atomic units, and currently assumes the number
  !! of atoms is the same as passed at construction/initialisation.
  !!
  subroutine send(this, energy, forces, stress)


    !> Instance
    class(IpiSocketComm), intent(inout) :: this


    !> Total energy
    real(dp), intent(in) :: energy


    !> Total forces
    real(dp), intent(in) :: forces(:,:)


    !> Cell stresses
    real(dp), intent(in) :: stress(3, 3)

    character(len=IPI_MSGLEN) :: header, buffer
    character(lc) :: msg
    integer :: nAtom

    @:ASSERT(this%tInit)
    @:ASSERT(size(forces, dim=1) == 3)

    nAtom = size(forces, dim=2)
    if (nAtom /= this%nAtom) then
      write(msg, '(1X,A,2I4)') 'Mismatch in number of atoms in socketSend', nAtom,this%nAtom
      call error(msg)
    end if

    call this%logger%write('ipisocket%send: Sending data to socket... ', 1)

    ! wait for anything other than a 'STATUS' state from the interface,
    ! returning state 'HAVEDATA' in the meanwhile
    listen: do while (.true.)
      call readbuffer(this%socket, header)
      call this%logger%write('ipisocket%send: read from socket: ' // trim(header), 3)

      if (trim(header)/='STATUS') then
        exit listen
      end if

      ! announce that we have available data
      buffer = 'HAVEDATA'
      call writebuffer(this%socket, buffer)
      call this%logger%write('ipisocket%send: write to socket: HAVEDATA', 3)
    end do listen

    ! expecting to send force data
    if (trim(header) /= 'GETFORCE') then
      call error("ipisocket%send: Error in socket communication! (expected 'GETFORCE',&
          & received '" // trim(header) // "'")
    end if

    buffer = 'FORCEREADY'
    call writebuffer(this%socket, buffer)
    call this%logger%write('ipisocket%send: write to socket: FORCEREADY', 3)

    ! transmit total energy
    call writebuffer(this%socket, energy)
    call this%logger%write('ipisocket%send: write to socket: energy', 3)
    write(msg, *) energy
    call this%logger%write(msg, 4)

    ! transmit number of atoms we have
    call writebuffer(this%socket, natom)
    call this%logger%write('ipisocket%send: write to socket: natom', 3)

    ! transmit forces
    call writebuffer(this%socket, reshape(forces, [3 * natom]))
    call this%logger%write('ipisocket%send: write to socket: forces', 3)
    call this%logger%write(forces, 5, '(f12.6)')

    ! transmit stress
    call writebuffer(this%socket, reshape(stress, [9]))
    call this%logger%write('ipisocket%send: write to socket: stress', 3)
    call this%logger%write(stress, 4, '(f12.6)')

    ! i-pi can also receive an arbitrary string, that will be printed out to the 'extra' trajectory
    ! file. this is useful if you want to return additional information, e.g.  atomic charges,
    ! wannier centres, etc. one must return the number of characters, then the string. here we just
    ! send back zero characters.
    call writebuffer(this%socket, 0)
    call this%logger%write('ipisocket%send: 0: nothing else to send', 3)
    call this%logger%write('ipisocket%send: Done', 1)

  end subroutine send


  !> Shuts down the socket.
  subroutine shutdown(this)


    !> Instance
    class(IpiSocketComm), intent(inout) :: this

    call close_socket(this%socket)
    this%tInit = .false.

  end subroutine shutdown

end module dftbp_ipisocket
