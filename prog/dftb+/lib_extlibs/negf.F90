!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"

!> Interfaces to libNEGF
module dftbp_extlibs_negf
  use libnegf, only : convertcurrent, eovh, getel, lnParams, pass_DM, Tnegf, kb, units
#:if WITH_MPI
  use libnegf, only : negf_mpi_init, negf_cart_init
#:endif
  use libnegf, only : z_CSR, z_DNS, READ_SGF, COMP_SGF, COMPSAVE_SGF
  use libnegf, only : DELTA_SQ, DELTA_W, DELTA_MINGO
  use libnegf, only : associate_lead_currents, associate_ldos, associate_transmission
  use libnegf, only : associate_current, compute_current, compute_density_dft, compute_ldos
  use libnegf, only : create, create_scratch, destroy, set_readoldDMsgf
  use libnegf, only : destroy_matrices, destroy_negf, get_params, init_contacts, init_ldos
  use libnegf, only : init_negf, init_structure, pass_hs, set_bp_dephasing
  use libnegf, only : set_drop, set_elph_block_dephasing, set_elph_dephasing, set_elph_s_dephasing
  use libnegf, only : set_ldos_indexes, set_params, set_scratch, set_tun_indexes, writememinfo
  use libnegf, only : writepeakinfo, dns2csr, csr2dns, nzdrop, printcsr
  use libnegf, only : compute_phonon_current, thermal_conductance
  use libnegf, only : convertHeatCurrent, convertHeatConductance
  implicit none

  public

end module dftbp_extlibs_negf
