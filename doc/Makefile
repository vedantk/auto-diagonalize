all: $(patsubst %.tex, %.pdf, $(wildcard *.tex))

%.pdf: %.tex
	pdflatex $^
