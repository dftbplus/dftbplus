!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"

#:assert not INSTANCE_SAFE_BUILD


!> Interfaces to libNEGF
module dftbp_extlibs_negf
  use libnegf, only : associate_current, associate_ldos, associate_lead_currents,&
      & associate_transmission, COMP_SGF, COMPSAVE_SGF, compute_current, compute_density_dft,&
      & compute_ldos, compute_phonon_current, convertcurrent, convertHeatConductance,&
      & convertHeatCurrent, create, create_scratch, csr2dns, DELTA_MINGO, DELTA_SQ, DELTA_W,&
      & destroy, destroy_matrices, destroy_negf, dns2csr, eovh, get_params, getel, init_contacts,&
      & init_ldos, init_negf, init_structure, kb, lnParams, nzdrop, pass_DM, pass_hs, printcsr,&
      & r_CSR, r_DNS, READ_SGF, set_bp_dephasing, set_drop, set_elph_block_dephasing,&
      & set_elph_dephasing, set_elph_s_dephasing, set_ldos_indexes, set_params, set_readoldDMsgf,&
      & set_scratch, set_tun_indexes, thermal_conductance, Tnegf, units, writememinfo,&
      & writepeakinfo, z_CSR, z_DNS
#:if WITH_MPI
  use libnegf, only : negf_cart_init, negf_mpi_init
#:endif
  implicit none

  public

end module dftbp_extlibs_negf
