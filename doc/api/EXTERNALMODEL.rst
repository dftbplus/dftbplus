Supplying external models to DFTB+
==================================

Use from within DFTB+::
  
  Model = External {
  
  }
  
Will use the externally declared model

Building
--------

Internally a toy model can be used, but if an external library is
provided in the configuration step (as a cmake argument to
statically/dynamically link a conforming external library) the
resulting libdftbplus or dftb+ can use the external model to generate
the hamiltonian or energy terms.

Model Units
-----------

All units are assumed to be atomic.


Declaring model capabilities
----------------------------

Initialising the model
----------------------

Updating the model
------------------

DFTB+ generates atom and bond-centred clusters surrounding the on-site
and diatomic parts of the model. Both types of cluster are surrounded
by `bystaner` atoms out to the `environment` cutoff (either around the
atom, or the bond vector between the atoms in the dimer). The dimers
of lengths of up to the `interaction` cutoff are generated.

In the special case of an environmentally independent model (i.e. the
environment cluster is of length 0), only individual atoms and pairs
of atoms will be generated.


Obtaining results
-----------------

Finalising and stopping
-----------------------
