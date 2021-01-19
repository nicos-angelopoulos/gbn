## gBN: BNs for large cohort genomic studies.

gBN runs the BN learning algorithm Gobnilp and implements many post-processing
visualisation functions for analysing large scale patient studies from genomic information.
Specially so for cancer patient data.

Currently only tested on *nix like systems. Developed on Linux Mint 2020.

## Executable dependencies

The following executables should be in the local path (see end of file for links).

* SWI-Prolog ()
* R
* gobnilp (SCIP-solver)
* scip ILP solver (needed by Gobnilp)

## R dependencies 

The following R packages will be asked for installation at loading time,
if not already present in local R installation.

* RColorBrewer
* cowplot
* ggplot2
* ggpubr
* gridExtra

## SWI pack dependencies

These will be installed (optionally) at loading time if they are not 
already present in local SWI installation.

* mtx
* real
* os\_lib
* by\_unix
* disp\_bn
* options
* debug\_call
* stoics\_lib
* svg

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


## pack info

* author nicos angelopoulos
* version  0.1 2020/9/26
* licence MIT

## Links
[gBN page](http://stoics.org.uk/~nicos/sware/gbn/)
[gBN git](https://github.com/nicos-angelopoulos/gbn)
[SWI-Prolog](https://www.swi-prolog.org/)
[R](https://www.r-project.org/)
[SCIP](https://scipopt.org/)
[GOBNILP](http://www.cs.york.ac.uk/aig/sw/gobnilp/)

## Author

Nicos Angelopoulos,\
London, 2020
