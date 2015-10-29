SHELL=/bin/bash  -i

v: all
	open faspsa.pdf

pdf:
	pdflatex faspsa.tex; bibtex faspsa; pdflatex faspsa.tex;

c clean clear:
	rm -f *.aux *.bbl *.blg *.log *.out *~
