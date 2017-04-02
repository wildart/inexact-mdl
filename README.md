# Inexact MDL for Linear Manifold Clusters

## Recreate figures

This LaTeX document contains embedded Julia code. In order to compile
the document correctly or extract source code, use Julia scientific report
generator package ``Weave.jl``.

### Requirements
- LaTex
  - pdflatex
- Julia 0.5
  - Packages: Weave, LMCLUS, DataFrames, Distributions, PLplot, Colors, Query, Gadfly

Use `jlweave` script to reproduce experiments and compile parts of the document, e.g.

	bin/jlweave results.texw

In order to recreate figures for the paper, run following commands:

    bin/jlweave src/results.texw
    bin/jlweave src/zero-dim-mdl.texw
    bin/jlweave src/mdl-error.texw
    bin/jlweave src/mdl-clust.texw
    julia src/climate_data_v3.jl


## Compile PDF

Run following command for a report PDF document

    bin/compile.sh report

Run following command to compile an ICPR paper PDF document:

    bin/compile.sh icpr

