## gBN: BNs for large cohort genomic studies.

gBN runs the BN learning algorithm Gobnilp and implements many post-processing
visualisation functions for analysing large scale patient studies from genomic information.
Specially so for cancer patient data.

Currently only tested on *nix like systems. Developed on Linux Mint 2020.\
It should also work on MacOS with one of the linux-like platforms installed on them.

## Executable dependencies

The following executables should be in the local path (see end of file for links).

* SWI-Prolog (swipl)
* R (compiled with --enable-R-shlib, the resulting libR.so should be in your LD_LIBRARY_PATH)
* gobnilp (C language version with SCIP-solver)
* scip ILP solver (needed by Gobnilp)
* graphviz (dot executable- for displaying)

## R dependencies 

The following R packages will be asked for installation at loading time,
if not already present in local R installation.

* RColorBrewer
* cowplot
* ggplot2
* ggpubr
* gridExtra

## SWI pack dependencies


Pack lib (pack(lib)), developed by us is central to loading all the Prolog dependencies.

```
% swipl
?- pack_install(lib).
```

The following will be installed, interactively, by pack(lib), at loading time if they are not 
already present in the local SWI installation. 

* mtx
* real
* os\_lib
* by\_unix
* disp\_bn
* options
* debug\_call
* stoics\_lib
* pack\_errors
* svg

You can install each pack above independently/separately from within SWI-Prolog, for example: 

```
% swipl
?- pack_install(mtx).
```


## Installation, Loading & Testing

```
% swipl

?- pack_install(gbn).

?- use_module(library(gbn)).
% Loading installed R library: RColorBrewer
% Loading installed R library: cowplot
% Loading installed R library: ggplot2
% Loading installed R library: ggpubr
% Loading installed R library: gridExtra

?- gbn([debug(true)]).
% Turning debugging on for predicate handle: gbn(gbn)
% Options: [$restore(gbn,debug,false),copy(false),data(pack(gbn/data/asia.dat)),display_dot(svg),odir(_33202),std_output(std_file)]
% Output directory: 'asia-20.09.27'
% Settings on: asia.set
```

A more complex example: 
```
% swipl

?- use_module(library(gbn)).

?- [cancer(aml)].

?- aml.
% Starting: aml
% Starting: fisher_nets
% Starting: fam_hmaps
% Starting: gates_nets
% Starting: svg_legend
% Finished: aml
true.
```

The above would create an output directory such as:  aml_min60-21.01.19/.

Cancer datasets from our paper (each can be run as per example above): 
* aml (Acute myeloid leukaemia)
* coa (Colon Adenocarcinoma, from TCGA)
* gbm (Glioblastoma)
* mpn (Myeloproliferative neoplasms)
* mye (Multiple myeloma)

## pack info

* author nicos angelopoulos
* version  0.1 2021/1/19
* license MIT

## Links
* [gBN page](http://stoics.org.uk/~nicos/sware/gbn/)
* [gBN git](https://github.com/nicos-angelopoulos/gbn)
* [SWI-Prolog](https://www.swi-prolog.org/)
* [R](https://www.r-project.org/)
* [SCIP](https://scipopt.org/)
* [GOBNILP](http://www.cs.york.ac.uk/aig/sw/gobnilp/)

## Author

Nicos Angelopoulos,\
London, 2021
