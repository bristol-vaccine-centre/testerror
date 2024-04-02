#!/bin/bash

cp refs.bib main/refs.bib
cp refs.bib s1/refs.bib
cp refs.bib s2/refs.bib

cd main
pdflatex  -interaction=batchmode main
bibtex main
pdflatex  -interaction=batchmode main
pdflatex  -interaction=batchmode main
latexdiff --math-markup=3 main.first_submission.tex main.tex > main.diff.tex
pdflatex  -interaction=batchmode main.diff
bibtex main.diff
pdflatex  -interaction=batchmode main.diff
pdflatex  -interaction=batchmode main.diff

cd ../s1
pdflatex  -interaction=batchmode s1
bibtex s1
pdflatex  -interaction=batchmode s1
pdflatex  -interaction=batchmode s1
latexdiff --math-markup=3 s1.first_submission.tex s1.tex > s1.diff.tex
pdflatex  -interaction=batchmode s1.diff
bibtex s1.diff
pdflatex  -interaction=batchmode s1.diff
pdflatex  -interaction=batchmode s1.diff

cd ../s2
pdflatex  -interaction=batchmode s2
bibtex s2
pdflatex  -interaction=batchmode s2
pdflatex  -interaction=batchmode s2
latexdiff --math-markup=3 s2.first_submission.tex s2.tex > s2.diff.tex
pdflatex  -interaction=batchmode s2.diff
bibtex s2.diff
pdflatex  -interaction=batchmode s2.diff
pdflatex  -interaction=batchmode s2.diff

