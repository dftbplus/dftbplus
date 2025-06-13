Types
=====

.. role:: f90(code)
   :language: Fortran

Enumerated cases
----------------

A number of modules use the following design pattern to provide a list
of the possible choices a variable can take 

.. code-block:: Fortran
   :linenos:

   !> Enumerator containing possible choices
   type :: TMyTypeEnum_

     !> Default case
     integer :: unknown = 0

     !> Case A
     integer :: A = 1

     !> Case two
     integer :: two = 2

   end type TMyTypeEnum_

   !> Actual instance of the MyType enumerator
   type(TMyType_), parameter :: myTypeList = TMyTypeEnum_()

Points to note are ::

  * The enumerated type is private and internal to the module (hence
    the underscore in the name), but the enumerated list is public, so
    visible in other modules.

  * The choces can be accessed as elements in the structure,
    :f90:`myTypeList%A`.

  * The variables of the enumerated type are all of the same kind.

Often this pattern is for a variable setting the state of an
associated type or class, called :f90:`type TMyType` in this example

.. code-block:: Fortran
   :linenos:

   !> Type for doing something
   type TMyType

     integer :: iMyType = MyTypeList%unknown

   end type TMyType

In this case, the initial state of the variable held in the type on
creation is `unknown`.
