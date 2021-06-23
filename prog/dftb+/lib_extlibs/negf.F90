!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"

#:assert not INSTANCE_SAFE_BUILD


!> Interfaces to libNEGF
module dftbp_extlibs_negf
  use libnegf, only : convertcurrent, eovh, getel, lnParams, pass_DM, Tnegf, kb, units,&
      & z_CSR, z_DNS, READ_SGF, COMP_SGF, COMPSAVE_SGF, DELTA_SQ, DELTA_W, DELTA_MINGO,&
      & associate_lead_currents, associate_ldos, associate_transmission, associate_current,&
      & compute_current, compute_density_dft, compute_ldos, create, create_scratch, destroy,&
      & set_readoldDMsgf, destroy_matrices, destroy_negf, get_params, init_contacts, init_ldos,&
      & init_negf, init_structure, pass_hs, set_bp_dephasing, set_drop, set_elph_block_dephasing,&
      & set_elph_dephasing, set_elph_s_dephasing, set_ldos_indexes, set_params, set_scratch,&
      & set_tun_indexes, writememinfo, writepeakinfo, dns2csr, csr2dns, nzdrop, printcsr,&
      & compute_phonon_current, thermal_conductance, convertHeatCurrent, convertHeatConductance
#:if WITH_MPI
  use libnegf, only : negf_mpi_init, negf_cart_init
#:endif
  implicit none

  public

end module dftbp_extlibs_negf
