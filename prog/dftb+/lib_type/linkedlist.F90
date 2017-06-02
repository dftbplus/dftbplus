!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains types and functions and subroutines for manipulating linked lists.
!!* Every list must be initialized with init, and destroyed with destroy.
module linkedList
#include "assert.h"
#include "allocate.h"  
  use accuracy
  implicit none
  private

  public :: listReal, listRealR1, listCharMc, listCharLc, listInt, listIntR1
  public :: init, destroy
  public :: append, len, find, hasElement, elemShape, isUnishaped
  public :: get, set, asArray, asVector, intoArray
  public :: charMc, charLc
  
  type listReal
    private
    integer                   :: length
    logical                   :: tUnishaped
    type(nodeReal),   pointer :: pFirst
    type(nodeReal),   pointer :: pLast
    integer                   :: iCache
    type(nodeReal),   pointer :: pCache
    logical                   :: tInitialized = .false.
  end type listReal

  type listRealR1
    private
    integer                   :: length
    integer                   :: elemShape(1)
    logical                   :: tUnishaped
    type(nodeRealR1), pointer :: pFirst
    type(nodeRealR1), pointer :: pLast
    integer                   :: iCache
    type(nodeRealR1), pointer :: pCache
    logical                   :: tInitialized = .false.
  end type listRealR1

  type listCharMc
    private
    integer                   :: length
    logical                   :: tUnishaped
    type(nodeCharMc), pointer :: pFirst
    type(nodeCharMc), pointer :: pLast
    integer                   :: iCache
    type(nodeCharMc), pointer :: pCache
    logical                   :: tInitialized = .false.
  end type listCharMc
  
  type listCharLc
    private
    integer                   :: length
    logical                   :: tUnishaped
    type(nodeCharLc), pointer :: pFirst
    type(nodeCharLc), pointer :: pLast
    integer                   :: iCache
    type(nodeCharLc), pointer :: pCache
    logical                   :: tInitialized = .false.
  end type listCharLc
  
  type listInt
    private
    integer                   :: length
    logical                   :: tUnishaped
    type(nodeInt),    pointer :: pFirst
    type(nodeInt),    pointer :: pLast
    integer                   :: iCache
    type(nodeInt),    pointer :: pCache
    logical                   :: tInitialized = .false.
  end type listInt

  type listIntR1
    private
    integer                   :: length
    integer                   :: elemShape(1)
    logical                   :: tUnishaped
    type(nodeIntR1),  pointer :: pFirst
    type(nodeIntR1), pointer  :: pLast
    integer                   :: iCache
    type(nodeIntR1), pointer  :: pCache
    logical                   :: tInitialized = .false.
  end type listIntR1

  !!* Generic interface for initializing lists
  interface init
    module procedure initReal
    module procedure initRealR1
    module procedure initCharMc
    module procedure initCharLc    
    module procedure initInt
    module procedure initIntR1
  end interface

  !!* Generic interface for destroying lists
  interface destroy
    module procedure destroyReal
    module procedure destroyRealR1
    module procedure destroyCharMc
    module procedure destroyCharLc
    module procedure destroyInt
    module procedure destroyIntR1
  end interface

  !!* Generic interface for appending elements to a list
  interface append
    module procedure appendReal
    module procedure appendRealR1
    module procedure appendCharMc
    module procedure appendCharLc
    module procedure appendInt
    module procedure appendIntR1
  end interface

  !!* Generic interface for getting the length(nr. of elements) of a list
  interface len
    module procedure lenReal
    module procedure lenRealR1
    module procedure lenCharMc
    module procedure lenCharLc
    module procedure lenInt
    module procedure lenIntR1
  end interface

  !!* Generic interface for getting the positon of an element in the list
  interface find
    module procedure findReal
    module procedure findRealR1
    module procedure findCharMc
    module procedure findCharLc
    module procedure findInt
    module procedure findIntR1
  end interface

  !!* Generic interface for checking if an element is in the list
  interface hasElement
    module procedure hasElementReal
    module procedure hasElementRealR1
    module procedure hasElementCharMc
    module procedure hasElementCharLc
    module procedure hasElementInt
    module procedure hasElementIntR1
  end interface
  
  !!* Generic interface for getting an element with a given index
  interface get
    module procedure getReal
    module procedure getRealR1
    module procedure getCharMc
    module procedure getCharLc
    module procedure getInt
    module procedure getIntR1
  end interface

  !!* Generic interface for getting an element with a given index
  interface set
    module procedure setReal
    module procedure setRealR1
    module procedure setCharMc
    module procedure setCharLc
    module procedure setInt
    module procedure setIntR1
  end interface

  !!* Generic interface for checking if all the list members have equal shape
  interface isUnishaped
    module procedure isUnishapedReal
    module procedure isUnishapedRealR1
    module procedure isUnishapedCharMc
    module procedure isUnishapedCharLc
    module procedure isUnishapedInt
    module procedure isUnishapedIntR1
  end interface
  
  !!* Generic interface for getting the list as an array
  interface asArray
    module procedure asArrayReal
    module procedure asArrayRealR1
    module procedure asArrayCharMc
    module procedure asArrayCharLc
    module procedure asArrayInt
    module procedure asArrayIntR1
  end interface

  !!* Generic interface for getting the list as a vector
  interface asVector
    module procedure asVectorRealR1
    module procedure asVectorIntR1
  end interface
  

  !!* Generic interface for getting the shape of the array
  interface elemShape
    module procedure getElemShapeRealR1
    module procedure getElemShapeIntR1
  end interface

  !! Generic interface to get rank one element in a bigger rank one array.
  interface intoArray
    module procedure intoArrayIntR1
    module procedure intoArrayRealR1
  end interface
  

  type nodeReal
    real(dp) :: value
    type(nodeReal), pointer :: pNext
  end type nodeReal

  type nodeRealR1
    real(dp),         pointer :: pValue(:)                    ! type specific
    type(nodeRealR1), pointer :: pNext
  end type nodeRealR1

  type nodeCharMc
    character(mc) :: value
    type(nodeCharMc), pointer :: pNext
  end type nodeCharMc
  
  type nodeCharLc
    character(lc) :: value
    type(nodeCharLc), pointer :: pNext
  end type nodeCharLc
  
  type nodeInt
    integer :: value
    type(nodeInt), pointer :: pNext
  end type nodeInt

  type nodeIntR1
    integer, pointer         :: pValue(:)
    type(nodeIntR1), pointer :: pNext
  end type nodeIntR1


contains

#include "real.inc"
#include "realr1.inc"
#include "charmc.inc"
#include "charlc.inc"  
#include "int.inc"
#include "intr1.inc"


  character(mc) function charMc(string)
    character(*), intent(in) :: string
    
    charMc = string(1:min(mc, len(string)))

  end function charMc

  

  character(lc) function charLc(string)
    character(*), intent(in) :: string
    
    charLc = string(1:min(lc, len(string)))

  end function charLc

  
end module linkedList
