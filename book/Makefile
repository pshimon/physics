NAME:=top
BLD_DIR=../build
$(NAME).pdf $(NAME).bcf: $(NAME).bbl $(NAME).tex
#	pdflatex $(NAME)
	xelatex --output-directory=$(BLD_DIR) $(NAME)
#	lualatex $(NAME)
$(NAME).bbl: $(NAME).bcf
	biber --output-directory=$(BLD_DIR)  $(NAME)

clean:
	rm -f $(NAME).pdf *.bbl *.bcf *.blg *.aux *.toc *.run.xml *.log *.nav *.out *.snm
