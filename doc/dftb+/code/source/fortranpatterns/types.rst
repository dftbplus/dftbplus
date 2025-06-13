Types
=====

.. role:: Fortran(code)
   :language: Fortran

Enumerated cases
----------------

A number of modules use the following design pattern to provide the
possible choices that a variable can take

.. code-block:: Fortran
   :linenos:

   !> Enumerator containing possible choices
   type :: TMyTypeEnum_

     !> Default case
     integer :: unknown = 0

     !> Case alpha
     integer :: alpha = 1

     !> Case beta
     integer :: beta = 2

   end type TMyTypeEnum_

   !> Actual instance of the MyType enumerator
   type(TMyType_), parameter :: myTypeEnum = TMyTypeEnum_()

Points to note ::

  * This type is private and internal to the module (hence the
    underscore in the name), but the resulting enumerated variable is
    public, so visible in other importing modules.

  * The available choices can be accessed as elements in the structure
    and their name should follow the camelCase convention.

  * The variables of the enumerated type are all of the same kind.

Often this pattern is for a variable setting the state of an
associated type or class, in this example the type would be called
:Fortran:`type TMyType` and the enumerated options referred to for
example as :Fortran:`if (ii /= myTypeEnum%unknown)`.

The associated type could look like

.. code-block:: Fortran
   :linenos:

   !> Type for doing something
   type TMyType

     integer :: iMyType = myTypeEnum%unknown

   end type TMyType

In this case, the initial state of the variable held in the type on
creation is `unknown`.
