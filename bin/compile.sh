#!/usr/bin/env bash

DOC=""
FINAL=0
BUILDDIR="build"
LATEXOPTS="-shell-escape -output-directory $BUILDDIR"
FILE=$1

if [ "$FILE" == "report" ]; then
    FINAL=1
    DOC="inexact_mdl_lmc"
elif [ "$FILE" == "icpr" ]; then
    FINAL=1
    DOC="icpr-2016"
elif [ "$FILE" = "draft" ]; then
    DOC="inexact_mdl_lmc"
elif [ "${FILE: -4}" == ".tex" ]; then
    DOC=`basename $FILE .tex`
    if [ ! -e $FILE ]; then
        echo "File $FILE does not exist"
        exit 1
    fi
else
    echo "$0 expected arguments: 'report', 'icpr', 'draft' or any tex file"
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
