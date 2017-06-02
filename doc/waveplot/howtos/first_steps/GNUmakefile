############################################################################
#
#  Simple tex file for creating howtos
#
############################################################################

# Adjust according the name of the howto
HOWTO = $(notdir $(CURDIR))

.SUFFIXES :
.SUFFIXES : .tex .pdf
.PHONY: all html pdf

HOWTO_PDF = $(patsubst %,%.pdf,$(HOWTO))
COMMONS = common
TEX_TEMP_SUFFIXES = .aux .out .log .toc
HTML_TEMP_SUFFIXES = .css
TEMPORARY_FILES = image.png labels.pl WARNINGS

PDFLATEX = /usr/bin/env pdflatex
LATEX2HTML = /usr/bin/env latex2html
LATEX2HTML_DESTDIR = html

all: pdf html

html: $(addsuffix .tex,$(HOWTO) $(COMMONS)) pdf
	$(LATEX2HTML) -init_file ../latex2html-init $<

pdf: $(HOWTO_PDF)

%.pdf: %.tex $(addsuffix .tex,$(COMMONS))
	$(PDFLATEX) $<
	$(PDFLATEX) $<


clean:
	rm -f $(notdir $(foreach suffix,$(TEX_TEMP_SUFFIXES),\
	  $(addsuffix $(suffix),$(HOWTO) $(COMMONS))))
	rm -f $(addprefix $(LATEX2HTML_DESTDIR)/, \
          $(foreach suffix,$(HTML_TEMP_SUFFIXES), \
          $(addsuffix $(suffix),$(HOWTO))))
	rm -f $(addprefix $(LATEX2HTML_DESTDIR)/,$(TEMPORARY_FILES))

distclean: clean
	rm -rf $(LATEX2HTML_DESTDIR)
	rm -f $(HOWTO_PDF)
