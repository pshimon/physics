NAME:=physics
BIBNAME:=phys-lib.bib
$(NAME).pdf $(NAME).bcf: $(NAME).bbl $(NAME).tex
	pdflatex $(NAME)
$(NAME).bbl: $(NAME).bcf
	biber $(NAME)

clean:
	rm -f $(NAME).pdf *.bbl *.bcf *.blg *.aux *.toc *.run.xml *.log *.nav *.out *.snm
