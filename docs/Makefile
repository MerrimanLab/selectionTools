PDFLATEX=pdflatex
LATEX_OPTS=-shell-escape
BIBTEXT=bibtex
INPUT_FILES=selection_pipeline.tex selection_pipeline.bib

all: $(INPUT_FILES)
	$(PDFLATEX) $(LATEX_OPTS) selection_pipeline.tex
	$(BIBTEXT) selection_pipeline
	$(PDFLATEX) $(LATEX_OPTS) selection_pipeline.tex
	$(PDFLATEX) $(LATEX_OPTS) selection_pipeline.tex	
clean:
	rm -rf selection_pipeline.pdf
