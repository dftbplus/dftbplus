#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2019  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

ROOT := $(PWD)

# Define a default goal here to make sure, make.config does not introduce on
.PHONY: _default
_default: default

include make.config

.PHONY: default misc all api

default: dftb+ modes waveplot
ifeq ($(strip $(WITH_TRANSPORT)),1)
  default: setupgeom
endif
misc: misc_skderivs misc_slakovalue
all: default misc
api: api_mm

.PHONY: install install_misc install_all install_api
install: install_dftb+ install_modes install_waveplot install_dptools
ifeq ($(strip $(WITH_TRANSPORT)),1)
  install: install_setupgeom
endif
install_misc: install_misc_skderivs install_misc_slakovalue install_tools_misc
install_all: install install_misc
install_api: install_api_mm

.PHONY: test test_api
test: test_dftb+ test_dptools

ifeq ($(strip $(BUILD_API)),1)
  test: test_api
endif
test_api: test_api_mm

.PHONY: check
check: check_dptools

################################################################################
# Sanity checks
################################################################################

# Check whether DEBUG level is correct
ifeq ($(filter 0 1 2,$(strip $(DEBUG))),)
  $(error 'Invalid value $(DEBUG) for DEBUG (must be 0, 1 or 2)')
endif

################################################################################
# Build targets
################################################################################

# You can disable automatic release determination by setting the release
# explicitely in the RELEASE file in the source root directory.
.PHONY: update_release
update_release:
	mkdir -p $(BUILDDIR)
	[ -r $(ROOT)/RELEASE ] && cp -a $(ROOT)/RELEASE $(BUILDDIR)/RELEASE \
        || $(ROOT)/utils/build/update_release $(BUILDDIR)/RELEASE \
        || echo "(UNKNOWN RELEASE)" > $(BUILDDIR)/RELEASE

.PHONY: dftb+ modes waveplot setupgeom
dftb+ modes waveplot setupgeom:
	mkdir -p $(BUILDDIR)/prog/$@
	$(MAKE) -C $(BUILDDIR)/prog/$@ -f $(ROOT)/prog/$@/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

DFTBPLUS_DEPS := update_release external_xmlf90
ifeq ($(strip $(WITH_SOCKETS)),1)
  DFTBPLUS_DEPS += external_fsockets
endif
ifeq ($(strip $(WITH_GPU)),1)
dftb+: external_magmahelper
endif
ifeq ($(strip $(WITH_DFTD3))$(strip $(COMPILE_DFTD3)),11)
  DFTBPLUS_DEPS += external_dftd3
endif
ifeq ($(strip $(WITH_MPI)),1)
  DFTBPLUS_DEPS += external_mpifx external_scalapackfx
endif
ifeq ($(strip $(WITH_TRANSPORT)),1)
  DFTBPLUS_DEPS += external_libnegf external_poisson
  external_poisson: external_libnegf
  ifeq ($(strip $(WITH_MPI)),1)
    external_libnegf: external_mpifx
    external_poisson: external_mpifx
  endif
endif
dftb+: $(DFTBPLUS_DEPS)

modes: external_xmlf90
waveplot: external_xmlf90
setupgeom: external_xmlf90

.PHONY: misc_skderivs misc_slakovalue
misc_skderivs misc_slakovalue:
	mkdir -p $(BUILDDIR)/prog/misc/$(subst misc_,,$@)
	$(MAKE) -C $(BUILDDIR)/prog/misc/$(subst misc_,,$@) \
	    -f $(ROOT)/prog/misc/$(subst misc_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

misc_skderivs misc_slakovalue: external_xmlf90


EXTERNAL_NAME = $(subst external_,,$@)

EXTERNALS = external_xmlf90 external_fsockets external_dftd3	\
    external_mpifx external_scalapackfx external_magmahelper	\
    external_poisson external_libnegf

.PHONY: $(EXTERNALS)
$(EXTERNALS):
	mkdir -p $(BUILDDIR)/external/$(EXTERNAL_NAME)
	$(MAKE) -C $(BUILDDIR)/external/$(EXTERNAL_NAME) \
          -f $(ROOT)/external/$(EXTERNAL_NAME)/make.dpbuild \
          ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)



.PHONY: api_mm
api_mm: api_lib_mm
ifeq ($(strip $(BUILD_TEST_BINARIES)),1)
  api_mm: api_tester_mm
endif


.PHONY: api_lib_mm
api_lib_mm:
	mkdir -p $(BUILDDIR)/api/mm
	$(MAKE) -C $(BUILDDIR)/api/mm -f $(ROOT)/api/mm/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

api_lib_mm: $(DFTBPLUS_DEPS)


API_TESTER_NAME = $(subst api_tester_,,$@)

.PHONY: api_tester_mm
api_tester_mm:
	mkdir -p $(BUILDDIR)/test/api/mm/testers
	$(MAKE) -C $(BUILDDIR)/test/api/mm/testers \
	    -f $(ROOT)/test/api/mm/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

api_tester_mm: api_lib_mm

################################################################################
# Test targets
################################################################################

.PHONY: test_dftb+
test_dftb+:
	$(MAKE) -C $(BUILDDIR)/prog/$(subst test_,,$@) \
	    -f $(ROOT)/prog/$(subst test_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR) test

test_dftb+: dftb+


.PHONY: test_dptools
test_dptools:
	mkdir -p $(BUILDDIR)/test/tools/dptools
	cd $< && $(ROOT)/test/tools/dptools/runtests.sh $(PYTHONS)


.PHONY: test_api_mm
test_api_mm:
	$(MAKE) -C $(BUILDDIR)/test/api/mm -f $(ROOT)/test/api/mm/make.build \
            ROOT=$(ROOT) BUILDROOT=$(BUILDDIR) test

test_api_mm: api_mm

################################################################################
# Install targets
################################################################################

.PHONY: install_dftb+ install_modes install_waveplot install_setupgeom
install_dftb+ install_modes install_waveplot install_setupgeom:
	$(MAKE) -C $(BUILDDIR)/prog/$(subst install_,,$@) \
	    -f $(ROOT)/prog/$(subst install_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR) install

install_dftb+: dftb+
install_modes: modes
install_waveplot: waveplot


.PHONY: install_misc_skderivs install_misc_slakovalue
install_misc_skderivs install_misc_slakovalue:
	$(MAKE) -C $(BUILDDIR)/prog/misc/$(subst install_misc_,,$@) \
	    -f $(ROOT)/prog/misc/$(subst install_misc_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR) install


PYTHON := python
.PHONY: install_dptools
install_dptools:
	cd $(ROOT)/tools/dptools \
            && $(PYTHON) setup.py install --prefix $(INSTALLDIR)


.PHONY: install_tools_misc
install_tools_misc:
	cp $(ROOT)/tools/misc/cubemanip $(INSTALLDIR)/bin
	cp $(ROOT)/tools/misc/plotxy $(INSTALLDIR)/bin


.PHONY: install_api_mm
install_api_mm:
	$(MAKE) -C $(BUILDDIR)/api/$(subst install_api_,,$@) \
	    -f $(ROOT)/api/$(subst install_api_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR) install

################################################################################
# Check targets
################################################################################
PYLINT2 := pylint
PYLINT3 := pylint3

.PHONY: check_dptools check_dptools_py2 check_dptools_py3
check_dptools: check_dptools_py2 check_dptools_py3
check_dptools_py2:
	$(PYLINT2) --rcfile utils/srccheck/pylint/pylintrc-2.ini\
            tools/dptools/src/dptools/ tools/dptools/bin/*[!~]
check_dptools_py3:
	$(PYLINT3) --rcfile utils/srccheck/pylint/pylintrc-3.ini\
            tools/dptools/src/dptools/ tools/dptools/bin/*[!~]

################################################################################
# Various targets
################################################################################

.PHONY: distclean
distclean:
	rm -rf $(BUILDDIR)

# Create a source distribution from current git check-out
# Note: check-out must contain all submodules
ARCHIVE_NAME := dftbplus
.PHONY: sourcedist
sourcedist:
	rm -rf $(BUILDDIR)/_sourcedist
	mkdir -p $(BUILDDIR)/_sourcedist
	$(ROOT)/utils/build/make_archive.sh $(ARCHIVE_NAME) \
            $(BUILDDIR)/_sourcedist
