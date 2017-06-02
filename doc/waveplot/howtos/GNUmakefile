############################################################################
#
#  Simple tex file for creating howtos
#
############################################################################
.SUFFIXES :
.SUFFIXES : .tex .html .pdf

HOWTOS = first_steps

TARGETS_CURDIR = package packageclean
TARGETS_HOWTODIR = $(filter-out $(TARGETS_CURDIR) $(HOWTOS),$(MAKECMDGOALS))

.PHONY: $(TARGETS_CURDIR) $(TARGETS_HOWTODIR) $(HOWTOS)

all: $(HOWTOS)

$(HOWTOS):
	@$(MAKE) -C $@

$(TARGETS_HOWTODIR):
	@for dir in $(HOWTOS); do $(MAKE) -C $$dir $@; done

PDFLATEX = /usr/bin/env pdflatex
LATEX2HTML = /usr/bin/env latex2html


.PHONY: package packageclean all
package: packageclean all
	mkdir -p $(addprefix package/,$(HOWTOS))
	for dir in $(HOWTOS); do \
          $(MAKE) -C $$dir clean; \
          cp -r $$dir/html package/$$dir; \
          if [ -e $$dir/download ]; then \
            cp -r $$dir/download package/$$dir/$$dir; \
            tar -C package/$$dir -c -z -f package/$$dir/$$dir.tar.gz $$dir; \
          fi; \
        done
	mkdir package/images
	cp images/*.png package/images
	cp howto.css package
	for howto in $(HOWTOS); do cp $$howto/$$howto.pdf package/$$howto; done

packageclean:
	rm -rf package
