% this has been transferred to prolog/gbn.pl, 28.08.16 SWI now displays all module predicates if that is defined prolog/gbn.pl

/** <module> Running and managing runs of gobnilp.


Run a GOBNILP, Bayesian networks (BNs) learning task and/or post-run routines on learning output.

Post processing includes the generation of Fisher based visualisation and family heatmaps.<br>
Assumes gobnlip is in your path. Currently only tested on *nix based systems.

---++ Installation

To install
==
?- pack_install(gbn).
==

to load
==
?- use_module(library(gbn)).
==

---++ Predicates

Main predicate for both running the BN experiment and post-processing is:
  * gbn/1

---++ Examples

==
?- gbn.
true.

?- ls.
% asia-24.07.13/   
true.
==

==
?- gbn(debug(true)).
% Turning debugging on for predicate handle: gbn(gbn)
% Options: [$restore(gbn,debug,false),copy(false),data(pack(gbn/data/asia.dat)),display_dot(svg),odir(_33202),std_output(std_file)]
% Output directory: 'asia-24.07.13'
% Settings on: asia.set
true.
==

A more complex example

==
?- [cancer(aml)].

?- absolute_file_name( pack('gbn/run/gbns_in_cancer'), Abs ),
|    ls( Abs ).
% aml.pl    coa.pl    data/     gbm.pl    mpn.pl    mye.pl    plots/    
Abs = '/home/nicos/.local/share/swi-prolog/pack/gbn/run/gbns_in_cancer'.

?- aml.
% Starting: aml
% Starting: fisher_nets
% Starting: fam_hmaps
% Starting: gates_nets
% Starting: svg_legend
% Finished: aml
true.
==

---++ Info
@author nicos angelopoulos
@version  0.0.1 2014/4/8
@version  0.1.0 2021/1/19
@version  0.2.0 2021/1/23
@see https://doi.org/10.1038/s42003-022-03243-w
@see https://stoics.org.uk/~nicos/sware/gbn
@see https://www.cs.york.ac.uk/aig/sw/gobnilp/
@see gbn/1
@see gbn_version/2

*/

/** gbn_module.

Documentation predicate to give anchor to the overall module documentation.

*/
gbn_module.
