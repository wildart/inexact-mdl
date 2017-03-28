Compile document
================

This LaTeX document contains embedded Julia code. In order to compile the document correctly or extract source code, use Julia scientific report generator ``Weave.jl``.

Required programs:
- LaTex
  - pdflatex, latexmk
- Julia
  - Packages: Weave, LMCLUS, DataFrames, Distributions, PLplot, Colors #Gadfly,

Use `jlweave` script to compile parts of the document:

	bin/jlweave results.texw

## Compile

To compile report run following commands

    bin/jlweave results.texw .
    bin/jlweave zero-dim-mdl.texw
    julia src/climate_data_v3.jl
    bin/compile.sh report

