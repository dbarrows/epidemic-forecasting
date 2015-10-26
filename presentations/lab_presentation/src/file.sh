#! /bin/bash
echo  > out.out
pdflatex slides
cp slides.pdf nonotes.pdf
echo \\setbeameroption{show only notes} > out.out