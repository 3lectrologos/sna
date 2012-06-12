.SUFFIXES: .pdf .dvi .tex
.PHONY: default clean distclean

REPORT = report

all: $(REPORT).pdf

%.pdf: %.tex *.bib
	latexmk -pdf $*
clean:
	@$(RM) *.aux *.lof *.log *.lot *.out *.toc *.nav *.gz *.pag *.snm *.bbl *.blg *~ \#*\# *.fdb_latexmk

distclean: clean
	@$(RM) $(REPORT).pdf $(REPORT).dvi $(REPORT).ps
