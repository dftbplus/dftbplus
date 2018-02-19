!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Storage of integer vectors in a FIFO (first in first out) queue and wrappers for different data
!> types and matrix ranks.
!> You can store and access vectors as FIFO. A given number of vectors is kept in the memory, if
!> there are more, they are written to disc.
!> note that in order to use the FIFO you have create and reset it.
module fifo
  use assert
  use accuracy
  use fileid
  implicit none
  private

  public :: OFifoIntR1, OFifoRealR1, OFifoRealR2, OFifoCplxR1, OFifoCplxR2
  public :: init, destruct, reset, get, push, restart


  !> First-in-first-out queue for reals
  type OFifoIntR1
    private

    !> Nr. of elements in stored vectors
    integer :: nElem

    !> Nr. of vectors kept in memory
    integer :: bufferSize

    !> File id for "swap" file
    integer :: fileId

    !> File name for "swap" file
    character(len=50) :: fileName

    !> Nr. of vectors stored in the fifo
    integer :: nStored

    !> Where to put next one in the buffer
    integer :: storePos

    !> Where to read the next from
    integer :: readPos

    !> Last operation (undefined/read/write)
    integer :: iMode

    !> Buffer full, open swap file
    logical :: tBufferFull

    !> The buffer itself.
    integer, allocatable :: buffer(:,:)

    !> Is the buffer initialised?
    logical :: tInit = .false.
  end type OFifoIntR1


  !> Real 1D vector
  type OFifoRealR1
    private
    type(OFifoIntR1) :: fifoIntR1
    integer :: nElem
    integer :: equivIntSize
    integer, allocatable :: convBuffer(:)
    logical :: tInit = .false.
  end type OFifoRealR1


  !> Real 2D array
  type OFifoRealR2
    private
    type(OFifoRealR1) :: fifoRealR1
    integer :: arrayShape(2)
    integer :: bufferSize
    logical :: tInit = .false.
  end type OFifoRealR2


  !> Complex 1D vector
  type OFifoCplxR1
    private
    type(OFifoIntR1) :: fifoIntR1
    integer :: nElem
    integer :: equivIntSize
    integer, allocatable :: convBuffer(:)
    logical :: tInit = .false.
  end type OFifoCplxR1


  !> Complex 2D array
  type OFifoCplxR2
    private
    type(OFifoCplxR1) :: fifoCplxR1
    integer :: arrayShape(2)
    integer :: bufferSize
    logical :: tInit = .false.
  end type OFifoCplxR2


  !> Initialises a fifo
  interface init
    module procedure FifoIntR1_init
    module procedure FifoRealR1_init
    module procedure FifoRealR2_init
    module procedure FifoCplxR1_init
    module procedure FifoCplxR2_init
  end interface init


  !> Destroys a fifo
  interface destruct
    module procedure FifoIntR1_destruct
    module procedure FifoRealR1_destruct
    module procedure FifoRealR2_destruct
    module procedure FifoCplxR1_destruct
    module procedure FifoCplxR2_destruct
  end interface destruct


  !> Resets the fifo
  interface reset
    module procedure FifoIntR1_reset
    module procedure FifoRealR1_reset
    module procedure FifoRealR2_reset
    module procedure FifoCplxR1_reset
    module procedure FifoCplxR2_reset
  end interface reset


  !> Returns the next vector from the FIFO.
  interface get
    module procedure FifoIntR1_get
    module procedure FifoRealR1_get
    module procedure FifoRealR2_get
    module procedure FifoCplxR1_get
    module procedure FifoCplxR2_get
  end interface get


  !> Puts a new vector on the FIFO.
  interface push
    module procedure FifoIntR1_push
    module procedure FifoRealR1_push
    module procedure FifoRealR2_push
    module procedure FifoCplxR1_push
    module procedure FifoCplxR2_push
  end interface push


  !> Rewinds for reading from the beginning
  interface restart
    module procedure FifoIntR1_restart
    module procedure FifoRealR1_restart
    module procedure FifoRealR2_restart
    module procedure FifoCplxR1_restart
    module procedure FifoCplxR2_restart
  end interface restart


  !> Undefined fifo state
  integer, parameter :: modeUndefined = 0

  !> Last operation was reading
  integer, parameter :: modeRead = 1

  !> Last operation was writing
  integer, parameter :: modeWrite = 2

contains


  !> Creates a FifoIntR1 instance, which stores rank 1 vectors.
  !> Caveat: this routine relies on the fact that file IDs have maximum of 5 digits.
  subroutine FifoIntR1_init(sf, bufferSize, fileName)

    !> Initialised FifoIntR1 instance on exit.
    type(OFifoIntR1), intent(out) :: sf

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in) :: bufferSize

    !> Name of the file to keep the vectors on the disc.
    character(len=*), intent(in) :: fileName

    @:ASSERT(bufferSize >= 0)
    @:ASSERT(len(fileName) > 0)

    sf%nElem = 0
    sf%bufferSize = bufferSize
    sf%fileId = getFileId()
    @:ASSERT(sf%fileId <= 99999)
    write(sf%fileName, "(A,I5.5)") &
        &fileName(1:min(len_trim(fileName),len(sf%fileName)-5)), sf%fileId
    allocate(sf%buffer(sf%nElem, sf%bufferSize))
    open(sf%fileId, file=sf%fileName, status="replace", form="unformatted",&
        &action="write")
    close(sf%fileId)
    sf%tInit = .true.

  end subroutine FifoIntR1_init


  !> Destroy FifoIntR1 object.
  subroutine FifoIntR1_destruct(sf)

    !> FifoIntR1 instance.
    type(OFifoIntR1), intent(inout) :: sf

    logical :: tOpened

    if (.not. sf%tInit) then
      return
    end if

    if (sf%tBufferFull) then
      inquire(sf%fileId, opened=tOpened)
      if (tOpened) then
        close(sf%fileId)
      end if
    end if
    open(sf%fileId, file=sf%fileName)
    close(sf%fileId, status="delete")

  end subroutine FifoIntR1_destruct


  !> Resets FifoIntR1
  subroutine FifoIntR1_reset(sf, nElem, bufferSize)

    !> FifoIntR1 instance.
    type(OFifoIntR1), intent(inout) :: sf

    !> Nr. of elements in the rank one vectors
    integer, intent(in) :: nElem

    !> In memory buffer size, before disc space is used
    integer, intent(in), optional :: bufferSize

    logical :: tOpened, tRealloc

    @:ASSERT(sf%tInit)
    @:ASSERT(nElem > 0)

    tRealloc = (nElem /= sf%nElem)
    if (present(bufferSize)) then
      if (bufferSize /= sf%bufferSize) then
        sf%bufferSize = bufferSize
        tRealloc = .true.
      end if
    end if

    ! Reallocate buffer if nr. of elements changed
    if (tRealloc) then
      sf%nElem = nElem
      deallocate(sf%buffer)
      allocate(sf%buffer(sf%nElem, sf%bufferSize))
    end if

    ! is there data left in the file on disc?
    if (sf%bufferSize < sf%nStored) then
      ! Empty file on the disc
      inquire(sf%fileId, opened=tOpened)
      if (tOpened) then
        close(sf%fileId)
      end if
      open(sf%fileId, file=sf%fileName, status="replace", form="unformatted",&
          &action="write")
      close(sf%fileId)
    end if

    ! Set initial variables
    if (sf%bufferSize == 0) then
      sf%tBufferFull = .true.
    else
      sf%tBufferFull = .false.
    end if
    sf%nStored = 0
    sf%storePos = 0
    sf%readPos = 0
    sf%iMode = modeUndefined

  end subroutine FifoIntR1_reset


  !> Adds a new vector to the FIFO.
  subroutine FifoIntR1_push(sf, vector)

    !> FifoIntR1 instance
    type(OFifoIntR1), intent(inout) :: sf

    !> Vector to store
    integer, intent(in) :: vector(:)

    integer :: ii
    logical :: tOpened

    @:ASSERT(sf%tInit)
    @:ASSERT(size(vector) == sf%nElem)

    ! Change to write mode if necessary
    if (sf%iMode /= modeWrite) then
      if (sf%tBufferFull) then
        inquire(sf%fileId, opened=tOpened)
        if (tOpened) then
          close(sf%fileId)
        end if
      end if
      sf%iMode = modeWrite
    end if

    sf%nStored = sf%nStored + 1

    ! If buffer size is zero, vector is written directly to the disc.
    if (sf%bufferSize == 0) then
      open(sf%fileId, file=sf%fileName, status="old", form="unformatted", &
          &position="append", action="write")
      write (sf%fileId) (vector(ii), ii = 1, sf%nElem)
      close(sf%fileId)
      return
    end if

    ! Get the next storing position /mod(storePos-1+1, bufferSize) + 1/
    sf%storePos = mod(sf%storePos, sf%bufferSize) + 1

    ! If buffer is full, oldest vector is written to disc before replaced
    if (sf%tBufferFull) then
      open(sf%fileId, file=sf%fileName, status="old", form="unformatted", &
          &position="append", action="write")
      write (sf%fileId) (sf%buffer(ii, sf%storePos), ii = 1, sf%nElem)
      close(sf%fileId)
    else
      sf%tBufferFull = (sf%nStored == sf%bufferSize)
    end if
    sf%buffer(:, sf%storePos) = vector(:)

  end subroutine FifoIntR1_push


  !> Returns the current vector from the FIFO without deleting it.
  !> This subroutines returns the vectors contained in the FIFO in the appropriate order (first in
  !> first out). The vectors are not deleted from the fifo. If all the vectors had been returned, or
  !> the subsequent get calls had been interrupted by a push-call, the vectors will be returned
  !> beginning with the most recent again.
  subroutine FifoIntR1_get(sf, vector)

    !> FifoIntR1 instance
    type(OFifoIntR1), intent(inout) :: sf

    !> Contains the vector on exit.
    integer, intent(out) :: vector(:)

    integer :: ind, ii

    @:ASSERT(sf%tInit)
    @:ASSERT(size(vector) == sf%nElem)
    @:ASSERT(sf%nStored > 0)

    if (sf%iMode /= modeRead) then
      call restart(sf)
    end if

    ! If vector to return is recent enough, return it from the buffer
    ! otherwise read it from the disc.
    sf%readPos = sf%readPos + 1
    if (sf%readPos > sf%nStored - sf%bufferSize) then
      ind = mod(sf%storePos + sf%bufferSize &
          &- (sf%nStored - sf%readPos) - 1, sf%bufferSize) + 1
      vector(:) = sf%buffer(:, ind)
    else
      read (sf%fileId) (vector(ii), ii = 1, sf%nElem)
    end if

    ! Close file if all old entries had been read from it.
    if (sf%tBufferFull &
        &.and. (sf%readPos == sf%nStored - sf%bufferSize)) then
      close(sf%fileId)
    end if

    ! If all entries had been returned, set fifo in an undefined state
    if (sf%readPos == sf%nStored) then
      sf%iMode = modeUndefined
    end if

  end subroutine FifoIntR1_get


  !> Rewinds the FIFO, so that reading starts with first element again.
  subroutine FifoIntR1_restart(sf)

    !> FIFO object
    type(OFifoIntR1), intent(inout) :: sf

    logical :: tOpened

    if (sf%tBufferFull) then
      inquire(sf%fileId, opened=tOpened)
      if (tOpened) then
        close(sf%fileId)
      end if
      open(sf%fileId, file=sf%fileName, status="old", form="unformatted",&
          &position="rewind", action="read")
    end if
    sf%readPos = 0
    sf%iMode = modeRead

  end subroutine FifoIntR1_restart


  !> Creates a FifoRealR1 instance, which stores rank 1 vectors.
  subroutine FifoRealR1_init(sf, bufferSize, fileName)

    !> Initialised FifoIntR1 instance on exit.
    type(OFifoRealR1), intent(out) :: sf

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in) :: bufferSize

    !> Name of the file to keep the vectors on the disc.
    character(len=*), intent(in) :: fileName

    @:ASSERT(.not. sf%tInit)

    call init(sf%fifoIntR1, bufferSize, fileName)
    sf%nElem = 0
    sf%equivIntSize = 0
    allocate(sf%convBuffer(sf%equivIntSize))
    sf%tInit = .true.

  end subroutine FifoRealR1_init


  !> Destruct FifoRealR1 object.
  subroutine FifoRealR1_destruct(sf)

    !> FifoRealR1 instance.
    type(OFifoRealR1), intent(inout) :: sf

    call destruct(sf%fifoIntR1)

  end subroutine FifoRealR1_destruct


  !> Resets FifoRealR1
  subroutine FifoRealR1_reset(sf, nElem, bufferSize)

    !> FifoRealR1 instance.
    type(OFifoRealR1), intent(inout) :: sf

    !> Nr. of elements in the rank one vectors
    integer, intent(in) :: nElem

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in), optional :: bufferSize

    real(dp), allocatable :: buffer(:)
    integer :: equiv(1)

    @:ASSERT(sf%tInit)
    @:ASSERT(nElem > 0)

    if (nElem /= sf%nElem) then
      deallocate(sf%convBuffer)
      sf%nElem = nElem
      allocate(buffer(nElem))
      sf%equivIntSize = size(transfer(buffer, equiv))
      deallocate(buffer)
      allocate(sf%convBuffer(sf%equivIntSize))
    end if
    if (present(bufferSize)) then
      call reset(sf%fifoIntR1, sf%equivIntSize, bufferSize)
    else
      call reset(sf%fifoIntR1, sf%equivIntSize)
    end if

  end subroutine FifoRealR1_reset


  !> Adds a new vector to the FIFO.
  subroutine FifoRealR1_push(sf, vector)

    !> FifoRealR1 instance
    type(OFifoRealR1), intent(inout) :: sf

    !> Vector to store
    real(dp), intent(in) :: vector(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(vector) == sf%nElem)

    sf%convBuffer = transfer(vector, sf%convBuffer, sf%equivIntSize)
    call push(sf%fifoIntR1, sf%convBuffer)

  end subroutine FifoRealR1_push


  !> Returns the current vector from the FIFO without deleting it.
  !> This subroutines returns the vectors contained in the FIFO in the appropriate order (first in
  !> first out). The vectors are not deleted from the fifo. If all the vectors had been returned, or
  !> the subsequent get calls had been interrupted by a push-call, the vectors will be returned
  !> beginning with the most recent again.
  subroutine FifoRealR1_get(sf, vector)

    !> FifoRealR1 instance
    type(OFifoRealR1), intent(inout) :: sf

    !> Contains the vector on exit.
    real(dp), intent(out) :: vector(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(vector) == sf%nElem)

    call get(sf%fifoIntR1, sf%convBuffer)
    vector(:) = transfer(sf%convBuffer, vector, sf%nElem)

  end subroutine FifoRealR1_get


  !> Rewinds the FIFO, so that reading starts with first element again.
  subroutine FifoRealR1_restart(sf)

    !> fifo object
    type(OFifoRealR1), intent(inout) :: sf

    call restart(sf%fifoIntR1)

  end subroutine FifoRealR1_restart


  !> Creates a FifoRealR2 instance, which stores rank 2 matrices.
  subroutine FifoRealR2_init(sf, bufferSize, fileName)

    !> Initialised FifoIntR2 instance on exit.
    type(OFifoRealR2), intent(out) :: sf

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in) :: bufferSize

    !> Name of the file to keep the vectors on the disc.
    character(len=*), intent(in) :: fileName

    @:ASSERT(.not. sf%tInit)

    sf%arrayShape(:) = 0
    sf%bufferSize = bufferSize
    call init(sf%fifoRealR1, 0, fileName)
    sf%tInit = .true.

  end subroutine FifoRealR2_init


  !> Destroy a FifoRealR2 object.
  subroutine FifoRealR2_destruct(sf)

    !> FifoRealR2 instance.
    type(OFifoRealR2), intent(inout) :: sf

    call destruct(sf%fifoRealR1)

  end subroutine FifoRealR2_destruct


  !> Resets FifoRealR2
  subroutine FifoRealR2_reset(sf, newShape)

    !> FifoRealR2 instance.
    type(OFifoRealR2), intent(inout) :: sf

    !> Nr. of elements in the vectors
    integer, intent(in) :: newShape(2)

    @:ASSERT(sf%tInit)
    @:ASSERT(all(newShape > 0))

    sf%arrayShape = newShape
    call reset(sf%fifoRealR1, sf%arrayShape(1), bufferSize=sf%bufferSize*sf%arrayShape(2))

  end subroutine FifoRealR2_reset


  !> Adds a new vector to the FIFO.
  subroutine FifoRealR2_push(sf, data)

    !> FifoRealR2 instance
    type(OFifoRealR2), intent(inout) :: sf

    !> Array to store
    real(dp), intent(in) :: data(:,:)

    integer :: ii

    @:ASSERT(sf%tInit)
    @:ASSERT(all(shape(data) == sf%arrayShape))

    do ii = 1, sf%arrayShape(2)
      call push(sf%fifoRealR1, data(:,ii))
    end do

  end subroutine FifoRealR2_push


  !> Returns the current vector from the FIFO without deleting it.
  !> This subroutines returns the vectors contained in the FIFO in the appropriate order (first in
  !> first out). The vectors are not deleted from the fifo. If all the vectors had been returned, or
  !> the subsequent get calls had been interrupted by a push-call, the vectors will be returned
  !> beginning with the most recent again.
  subroutine FifoRealR2_get(sf, data)

    !> FifoRealR2 instance
    type(OFifoRealR2), intent(inout) :: sf

    !> Contains the vector on exit.
    real(dp), intent(out) :: data(:,:)

    integer :: ii

    @:ASSERT(sf%tInit)
    @:ASSERT(all(shape(data) == sf%arrayShape))

    do ii = 1, sf%arrayShape(2)
      call get(sf%fifoRealR1, data(:,ii))
    end do

  end subroutine FifoRealR2_get


  !> Rewinds the FIFO, so that reading starts with first element again.
  subroutine FifoRealR2_restart(sf)

    !> FifoRealR2 instance
    type(OFifoRealR2), intent(inout) :: sf

    call restart(sf%fifoRealR1)

  end subroutine FifoRealR2_restart


  !> Creates a FifoCplxR1 instance, which stores rank 1 vectors.
  subroutine FifoCplxR1_init(sf, bufferSize, fileName)

    !> Initialised FifoIntR1 instance on exit.
    type(OFifoCplxR1), intent(out) :: sf

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in) :: bufferSize

    !> Name of the file to keep the vectors on the disc.
    character(len=*), intent(in) :: fileName

    @:ASSERT(.not. sf%tInit)

    call init(sf%fifoIntR1, bufferSize, fileName)
    sf%nElem = 0
    sf%equivIntSize = 0
    allocate(sf%convBuffer(sf%equivIntSize))
    sf%tInit = .true.

  end subroutine FifoCplxR1_init


  !> Destroy a FifoCplxR1 object.
  subroutine FifoCplxR1_destruct(sf)

    !> FifoCplxR1 instance.
    type(OFifoCplxR1), intent(inout) :: sf

    call destruct(sf%fifoIntR1)

  end subroutine FifoCplxR1_destruct


  !> Resets a FifoCplxR1
  subroutine FifoCplxR1_reset(sf, nElem, bufferSize)

    !> FifoCplxR1 instance.
    type(OFifoCplxR1), intent(inout) :: sf

    !> Nr. of elements in the rank one vectors
    integer, intent(in) :: nElem

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in), optional :: bufferSize

    complex(dp), allocatable :: buffer(:)
    integer :: equiv(1)

    @:ASSERT(sf%tInit)
    @:ASSERT(nElem > 0)

    if (nElem /= sf%nElem) then
      deallocate(sf%convBuffer)
      sf%nElem = nElem
      allocate(buffer(nElem))
      sf%equivIntSize = size(transfer(buffer, equiv))
      deallocate(buffer)
      allocate(sf%convBuffer(sf%equivIntSize))
    end if
    if (present(bufferSize)) then
      call reset(sf%fifoIntR1, sf%equivIntSize, bufferSize)
    else
      call reset(sf%fifoIntR1, sf%equivIntSize)
    end if

  end subroutine FifoCplxR1_reset


  !> Adds a new vector to the FIFO.
  subroutine FifoCplxR1_push(sf, vector)

    !> FifoCplxR1 instance
    type(OFifoCplxR1), intent(inout) :: sf

    !> Vector to store
    complex(dp), intent(in) :: vector(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(vector) == sf%nElem)

    sf%convBuffer = transfer(vector, sf%convBuffer, sf%equivIntSize)
    call push(sf%fifoIntR1, sf%convBuffer)

  end subroutine FifoCplxR1_push


  !> Returns the current vector from the FIFO without deleting it.
  !> This subroutines returns the vectors contained in the FIFO in the appropriate order (first in
  !> first out). The vectors are not deleted from the fifo. If all the vectors had been returned, or
  !> the subsequent get calls had been interrupted by a push-call, the vectors will be returned
  !> beginning with the most recent again.
  subroutine FifoCplxR1_get(sf, vector)

    !> FifoCplxR1 instance
    type(OFifoCplxR1), intent(inout) :: sf

    !> Contains the vector on exit.
    complex(dp), intent(out) :: vector(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(vector) == sf%nElem)

    call get(sf%fifoIntR1, sf%convBuffer)
    vector(:) = transfer(sf%convBuffer, vector, sf%nElem)

  end subroutine FifoCplxR1_get


  !> Rewinds the FIFO, so that reading starts with first element again.
  subroutine FifoCplxR1_restart(sf)

    !> FifoCplxR1 instance
    type(OFifoCplxR1), intent(inout) :: sf

    call restart(sf%fifoIntR1)

  end subroutine FifoCplxR1_restart


  !> Creates a FifoCplxR2 instance, which stores rank 2 matrices.
  subroutine FifoCplxR2_init(sf, bufferSize, fileName)

    !> Initialised FifoCplxR2 instance on exit.
    type(OFifoCplxR2), intent(out) :: sf

    !> Nr. of vectors to keep in the memory (>=0)
    integer, intent(in) :: bufferSize

    !> Name of the file to keep the vectors on the disc.
    character(len=*), intent(in) :: fileName

    @:ASSERT(.not. sf%tInit)

    sf%arrayShape(:) = 0
    sf%bufferSize = bufferSize
    call init(sf%fifoCplxR1, 0, fileName)
    sf%tInit = .true.

  end subroutine FifoCplxR2_init


  !> Destroy a FifoCplxR2 object.
  subroutine FifoCplxR2_destruct(sf)

    !> FifoCplxR2 instance.
    type(OFifoCplxR2), intent(inout) :: sf

    call destruct(sf%fifoCplxR1)

  end subroutine FifoCplxR2_destruct


  !> Resets a FifoCplxR2 instance
  subroutine FifoCplxR2_reset(sf, newShape)

    !> FifoCplxR2 instance.
    type(OFifoCplxR2), intent(inout) :: sf

    !> Nr. of elements in the arrays
    integer, intent(in) :: newShape(2)

    @:ASSERT(sf%tInit)
    @:ASSERT(all(newShape > 0))

    sf%arrayShape = newShape
    call reset(sf%fifoCplxR1, sf%arrayShape(1), bufferSize=sf%bufferSize*sf%arrayShape(2))

  end subroutine FifoCplxR2_reset


  !> Adds a new vector to the FIFO.
  subroutine FifoCplxR2_push(sf, data)

    !> FifoCplxR2 instance
    type(OFifoCplxR2), intent(inout) :: sf

    !> Vector to store
    complex(dp), intent(in) :: data(:,:)

    integer :: ii

    @:ASSERT(sf%tInit)
    @:ASSERT(all(shape(data) == sf%arrayShape))

    do ii = 1, sf%arrayShape(2)
      call push(sf%fifoCplxR1, data(:,ii))
    end do

  end subroutine FifoCplxR2_push


  !> Returns the current vector from the FIFO without deleting it.
  !> This subroutines returns the vectors contained in the FIFO in the appropriate order (first in
  !> first out). The vectors are not deleted from the fifo. If all the vectors had been returned, or
  !> the subsequent get calls had been interrupted by a push-call, the vectors will be returned
  !> beginning with the most recent again.
  subroutine FifoCplxR2_get(sf, data)

    !> FifoCplxR2 instance
    type(OFifoCplxR2), intent(inout) :: sf

    !> Contains the array on exit.
    complex(dp), intent(out) :: data(:,:)

    integer :: ii

    @:ASSERT(sf%tInit)
    @:ASSERT(all(shape(data) == sf%arrayShape))

    do ii = 1, sf%arrayShape(2)
      call get(sf%fifoCplxR1, data(:,ii))
    end do

  end subroutine FifoCplxR2_get


  !> Rewinds the FIFO, so that reading starts with first element again.
  subroutine FifoCplxR2_restart(sf)

    !> instance of a fifo
    type(OFifoCplxR2), intent(inout) :: sf

    call restart(sf%fifoCplxR1)

  end subroutine FifoCplxR2_restart

end module fifo
