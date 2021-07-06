## gBN: BNs for large cohort genomic studies.

gBN runs the BN learning algorithm Gobnilp and implements many post-processing
visualisation functions for analysing large scale patient studies from genomic information.
Specially so for cancer patient data.

Currently only tested on *nix like systems. Developed on Linux Mint 2020.\
It should also work on MacOS with one of the linux-like platforms installed.

## Executable dependencies

The following executables should be installed on your local machine.

* [SWI-Prolog](https://www.swi-prolog.org/) (swipl in $PATH)
* [R](https://www.r-project.org/) (compiled with --enable-R-shlib, the resulting _libR.so_ should be in $LD_LIBRARY_PATH)
* [GOBNILP](http://www.cs.york.ac.uk/aig/sw/gobnilp/) (C language version with SCIP-solver; _gobnilp_ should be in $PATH)
* [SCIP](https://scipopt.org/) SCIP ILP solver (needed by Gobnilp at installation time)
* [graphviz](https://graphviz.org/) (for displaying, _dot_ executable should be in $PATH)

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

All packs installable this way can be found on the SWI-Prolog [pack list](https://www.swi-prolog.org/pack/list).

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

The above would create an output directory: aml_min60-21.01.19 where the date stamp will reflect
the current date.

Cancer datasets from our paper (each can be run as per example above): 
* aml (Acute myeloid leukaemia)
* coa (Colon Adenocarcinoma, from TCGA)
* gbm (Glioblastoma)
* mpn (Myeloproliferative neoplasms)
* mye (Multiple myeloma)

## Data

Datasets are in the source directory: data/ (with cancer datasets in data/gbns_in_cancer/).

## Documentation

* [static lib documentation](http://stoics.org.uk/~nicos/sware/gbn/doc/html/gbn.html)

## Main predicates - worked example

The main flow of execution is to create a run of Gobnilp experiments which are deposited in a new directory.
These create the vanilla Gobnilp outputs. gbn predicates are then called to operate on this directory to 
create the additional analyses and outputs.

For example in a new empty directory start an SWI session and load the library (here shown in Linux):
```
> mkdir /tmp/gbn_example
> cd /tmp/gbn_example
> swipl 

Welcome to SWI-Prolog (threaded, 64 bits, version 8.3.26-20-g3cd115ec2)
SWI-Prolog comes with ABSOLUTELY NO WARRANTY. This is free software.
Please run ?- license. for legal details.

For online help and background, visit https://www.swi-prolog.org
For built-in help, use ?- help(Topic). or ?- apropos(Word).

?- use_module(library(gbn)).
% Loading installed R library: RColorBrewer
% Loading installed R library: cowplot
% Loading installed R library: ggplot2
% Loading installed R library: ggpubr
% Loading installed R library: gridExtra
true.

?- 
```

The first step in a gbn analysis is to run Gobnilp on some data with predicate gbn/1. This takes a list as argument, 
the items of which: (a) point to the dataset is to be analysed, (b) where the results are to be placed and (c) change
Gobnilp settings for this run.

```
?- GBNOpts = [ data(pack('gbn/data/gbns_in_cancer/aml/aml_min60.dat')),
          dir(OsOdir),
          setting(edge_penalty,7)
        ],
   gbn( GBNOpts ).

GBNOpts = [data(pack('gbn/data/gbns_in_cancer/aml/aml_min60.dat')), dir('aml_min60-21.07.06'), setting(edge_penalty, 7)],
OsOdir = 'aml_min60-21.07.06'.

```
In the above example the aml_min_60.dat from within the gbn installation is used.
Output directory is given a variable, which means gbn creates a default name. In this case: OsOdir = 'aml_min60-21.07.06'.
The option setting(edge_penalty,7) instructs Gobnilp to use value 7 for its edge_penalty setting. 

Note, that often a directory will hold a number of Gobnilp runs on a single dataset where a number of settings might be tried.
For the instance the following call will run 10 Gobnilp runs, each with a distinct integer value for edge_penalty in the
range of 1-10.

```
?- GBNOpts = [ data(pack('gbn/data/gbns_in_cancer/aml/aml_min60.dat')),
          dir(OsOdir),
          multiple(edge_penalty,[1,2,3,4,5,6,7,8,9,10])
        ],
   gbn( GBNOpts ).

GBNOpts = [data(pack('gbn/data/gbns_in_cancer/aml/aml_min60.dat')), dir('aml_min60-21.07.06.12.31'), multiple(edge_penalty, [1, 2, 3, 4, 5|...])],
OsOdir = 'aml_min60-21.07.06.12.31'.
```

The new directory contents are:
```
?- ls('aml_min60-21.07.06.12.31').
% aml_min60.dat                 aml_min60-e2.bn               aml_min60-e4.dot              aml_min60-e6.set              aml_min60-e8.svg
% aml_min60-e10.bn              aml_min60-e2.dot              aml_min60-e4.set              aml_min60-e6.svg              aml_min60-e9.bn
% aml_min60-e10.dot             aml_min60-e2.set              aml_min60-e4.svg              aml_min60-e7.bn               aml_min60-e9.dot
% aml_min60-e10.set             aml_min60-e2.svg              aml_min60-e5.bn               aml_min60-e7.dot              aml_min60-e9.set
% aml_min60-e10.svg             aml_min60-e3.bn               aml_min60-e5.dot              aml_min60-e7.set              aml_min60-e9.svg
% aml_min60-e1.bn               aml_min60-e3.dot              aml_min60-e5.set              aml_min60-e7.svg              aml_min60_scnti.txt
% aml_min60-e1.dot              aml_min60-e3.set              aml_min60-e5.svg              aml_min60-e8.bn               aml_min60.set
% aml_min60-e1.set              aml_min60-e3.svg              aml_min60-e6.bn               aml_min60-e8.dot              aml_min60_std_output.txt
% aml_min60-e1.svg              aml_min60-e4.bn               aml_min60-e6.dot              aml_min60-e8.set              
true.

?- 

```

Given a directory of Gobnilp experiments on a single dataset, pack(gbn) has a number of predicates that can post process this directory to produce a number 
of additional outputs.

gbn_fisher_nets/1 can be used to produce visualisations of the produced BNs where edges are coloured according to Fisher statistic, where colour depends 
on mutual-exclusivity or co-occurance and continuity of the line shows significance or not (doted line). In the single edge penalty example above 
(edge_penalty=7) the following shos the contents of the output directory before and after the execution of gbn_fisher_nets/1.

```
?- ls('aml_min60-21.07.06').
% aml_min60.bn                aml_min60.dot               aml_min60.set               aml_min60.svg               
% aml_min60.dat               aml_min60_scnti.txt         aml_min60_std_output.txt    
true.

?- gbn_fisher_nets(dir('aml_min60-21.07.06')).
true.

?- ls('aml_min60-21.07.06').
% aml_min60.bn                aml_min60.dot               aml_min60_fclr.dot          aml_min60_fisher.dot        aml_min60_scnti.txt         aml_min60_std_output.txt
% aml_min60.dat               aml_min60_fclr.csv          aml_min60_fclr.svg          aml_min60_fisher.svg        aml_min60.set               aml_min60.svg
true.

?- 
```

gbn_fam_hmaps/1 produces a heatmap for each family in each output BN within the output directory is given to work on. In addition, the predicate generates single 
pdf collecting all family heatmaps. These can be extremely useful for including into supplementary information of publications. Continue on the above example, we have:

```
?- gbn_fam_hmaps(dir('aml_min60-21.07.06')).
true.

?- ls('aml_min60-21.07.06').
% aml_min60.bn                aml_min60_fams/             aml_min60_fclr.svg          aml_min60_gates_best.csv    aml_min60_scnti.txt         aml_min60.svg
% aml_min60.dat               aml_min60_fclr.csv          aml_min60_fisher.dot        aml_min60_multi_prns.png    aml_min60.set               
% aml_min60.dot               aml_min60_fclr.dot          aml_min60_fisher.svg        aml_min60_prns/             aml_min60_std_output.txt    
true.
```

## Raspberry pi 4 OS image

We also provide a full Raspberry pi 4 OS which has all dependencies installed within an Ubuntu 20.10 operating system.

* [gBN Rasp4 image](http://stoics.org.uk/~nicos/sware/gbn/gbn_image.html)

The compressed image is 3Gb in size and it will uncompress to 14Gb in size.

## Pack info

* author nicos angelopoulos
* version  0.1 2021/1/19
* license MIT
* [gBN page](http://stoics.org.uk/~nicos/sware/gbn/)
* [gBN git](https://github.com/nicos-angelopoulos/gbn)
* [gBN at SWI-Prolog](https://www.swi-prolog.org/pack/list?p=gbn)

## Publications

* Classification and Personalized Prognosis in Myeloproliferative Neoplasms, [10.1056/NEJMoa1716614](http://dx.doi.org/10.1056/NEJMoa1716614)
* Genomic Landscape and Chronological Reconstruction of Driver Events in Multiple Myeloma, [10.1038/s41467-019-11680-1](https://doi.org/10.1038/s41467-019-11680-1)
* Molecular Evolution of IDH Wild-Type Glioblastomas Treated With Standard of Care Affects Survival and Design of Precision Medicine Trials: A Report From the EORTC 1542 Study, [10.1200/JCO.19.00367](https://ascopubs.org/doi/abs/10.1200/JCO.19.00367)

## Author

Nicos Angelopoulos,\
London, 2021
