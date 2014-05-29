# Make file for making manual
#
# author James Boocock

PDFLATEX=pdflatex
LATEX_OPTS=-shell-escape
BIBTEXT=bibtex
INPUT_FILES=selection_pipeline.tex selection_pipeline.bib

all: $(INPUT_FILES)
	$(PDFLATEX) $(LATEX_OPTS) selection_pipeline.tex
	$(BIBTEXT) selection_pipeline
	$(PDFLATEX) $(LATEX_OPTS) selection_pipeline.tex
	$(PDFLATEX) $(LATEX_OPTS) selection_pipeline.tex	
	rm -f selection_pipeline.aux
	rm -f selection_pipeline.bbl
	rm -f selection_pipeline.blg
	rm -f selection_pipeline.dvi
	rm -f selection_pipeline.idx
	rm -f selection_pipeline.toc
	rm -f selection_pipeline.out
	rm -f selection_pipeline.log
clean:
	rm -rf selection_pipeline.pdf
