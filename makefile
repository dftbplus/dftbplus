#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

ROOT := $(PWD)

.PHONY: default misc all
default: dftb+ modes waveplot
misc: misc_skderivs misc_slakovalue
all: default misc

.PHONY: install install_misc install_all
install: install_dftb+ install_modes install_waveplot install_dptools
install_misc: install_misc_skderivs install_misc_slakovalue

.PHONY: test
test: test_dftb+ test_dptools

.PHONY: check
check: check_dptools

include make.config

################################################################################
# Sanity checks
################################################################################

# Check whether DEBUG level is correct
ifeq ($(filter 0 1 2,$(strip $(DEBUG))),)
  $(error 'Invalid value $(DEBUG) for DEBUG (must be 0, 1 or 2)')
endif

# Check whether PROGRESS use serial compilation
ifeq ($(strip $(WITH_PROGRESS))$(strip $(WITH_MPI)),11)
  $(error 'The PROGRESS library can not be used in MPI binaries')
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

.PHONY: dftb+ modes waveplot
dftb+ modes waveplot:
	mkdir -p $(BUILDDIR)/prog/$@
	$(MAKE) -C $(BUILDDIR)/prog/$@ -f $(ROOT)/prog/$@/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

dftb+: update_release external_xmlf90
ifeq ($(strip $(WITH_SOCKETS)),1)
dftb+: external_fsockets
endif
ifeq ($(strip $(WITH_DFTD3))$(strip $(COMPILE_DFTD3)),11)
dftb+: external_dftd3
endif
ifeq ($(strip $(WITH_MPI)),1)
dftb+: external_mpifx external_scalapackfx
endif
modes: external_xmlf90
waveplot: external_xmlf90


.PHONY: misc_skderivs misc_slakovalue
misc_skderivs misc_slakovalue:
	mkdir -p $(BUILDDIR)/prog/misc/$(subst misc_,,$@)
	$(MAKE) -C $(BUILDDIR)/prog/misc/$(subst misc_,,$@) \
	    -f $(ROOT)/prog/misc/$(subst misc_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

misc_skderivs: external_xmlf90


EXTERNAL_NAME = $(subst external_,,$@)

EXTERNALS = external_xmlf90 external_fsockets external_dftd3 external_mpifx\
    external_scalapackfx
.PHONY: $(EXTERNALS)
$(EXTERNALS):
	mkdir -p $(BUILDDIR)/external/$(EXTERNAL_NAME)
	$(MAKE) -C $(BUILDDIR)/external/$(EXTERNAL_NAME) \
          -f $(ROOT)/external/$(EXTERNAL_NAME)/make.dpbuild \
          ROOT=$(ROOT) BUILDROOT=$(BUILDDIR)

################################################################################
# Test targets
################################################################################

.PHONY: test_dftb+
test_dftb+:
	$(MAKE) -C $(BUILDDIR)/prog/$(subst test_,,$@) \
	    -f $(ROOT)/prog/$(subst test_,,$@)/make.build \
	    ROOT=$(ROOT) BUILDROOT=$(BUILDDIR) test

test_dftb+: dftb+


test_dptools:
	mkdir -p $(BUILDDIR)/test/tools/dptools
	cd $< && $(ROOT)/test/tools/dptools/runtests.sh $(PYTHONS)


################################################################################
# Install targets
################################################################################

.PHONY: install_dftb+ install_modes install_waveplot
install_dftb+ install_modes install_waveplot:
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
