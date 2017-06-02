.SUFFIXES :
.SUFFIXES : .fig .pdf .ps .toc .bbl .aux .tex .ind .html ._tmp_

MANUALCOMPONENTS = manual.tex

RM = /usr/bin/env rm -f
FIG2DEV = /usr/bin/env fig2dev
PDFOPTSFIG2DEV = -L pdf
PS2PDF = /usr/bin/env ps2pdf
PDF2PS = /usr/bin/env pdf2ps
PDF2PSOPTS = -sPAPERSIZE=a4
LATEX = /usr/bin/env pdflatex
BIBTEX = /usr/bin/env bibtex
MAKEINDEX = /usr/bin/env makeindex
LATEX2HTML = /usr/bin/env latex2html

# default rule
default: manual

%.pdf: %.fig  #convert .fig files into .pdf
	${FIG2DEV} ${PDFOPTSFIG2DEV} $< $*.pdf

%.pdf: %.tex  # TEX to PDF
	${LATEX} $*
	${LATEX} $*

.PRECIOUS: %._tmp_
%._tmp_: %.tex
	${LATEX} $*
	touch $@

%.toc: %.tex
	${LATEX} $*

%.ps: %.pdf    # PDF 2 PS
	${PDF2PS} $< $@

%.bbl: %._tmp_ %.bib
	${BIBTEX} $*

%.ind: %._tmp_
	${MAKEINDEX} $*

%.html: %._tmp_
	${LATEX2HTML} $*

manual : manual.pdf

manual._tmp_: ${MANUALCOMPONENTS}

manual.pdf : manual.toc manual.ind

manual.ps : manual.pdf

manual.html : manual.tex

dependancies.mk : manual.tex
	./depend.pl manual.tex > dependancies.mk

clean :
	${RM} manual.pdf manual.toc manual.out *.aux \
  manual.log manual.blg manual.brf manual.ps dependancies.mk \
  manual.ind manual.idx manual.ilg manual._tmp_
	${RM} -r manual

