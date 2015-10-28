SHELL=/bin/bash  -i

v: all 
	open report.pdf

all:
	pdflatex report.tex; bibtex report; pdflatex report.tex;

clear:
	rm -f *.aux *.bbl *.blg *.log *.out *~

c:
	rm -f *.aux *.bbl *.blg *.log *.out *.pdf *~

