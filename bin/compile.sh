#!/bin/sh

DOC="inexact_mdl_lmc"
FINAL=1
BUILDDIR="build"
LATEXOPTS="-shell-escape -output-directory $BUILDDIR"

if [ "$1" = "report" ]; then
    DOC="inexact_mdl_lmc"
elif [ "$1" = "icpr" ]; then
    DOC="icpr-2016"
elif [ "$1" = "draft" ]; then
    FINAL=0
else
    echo "$0 expected arguments: 'report', 'icpr' or 'draft'"
    exit 1
fi

mkdir -p $BUILDDIR

pdflatex $LATEXOPTS $DOC.tex
if [ $FINAL = 1 ]; then
    bibtex $BUILDDIR/$DOC.aux
    pdflatex $LATEXOPTS $DOC.tex
    pdflatex $LATEXOPTS $DOC.tex
fi

cp $BUILDDIR/$DOC.pdf .

# Cleanup
if [ $FINAL = 1 ]; then
    rm -rf build
fi
# for EXT in "aux" "log" "out" "toc" "bbl" "blg"
# do
#     rm $DOC.$EXT
# done
