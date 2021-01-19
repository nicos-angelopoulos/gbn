/** <module> Running and managing runs of gobnilp.

	Run the BN learning algorithm gobnilp.

	Assumes gobnlip is in your path. Currently only tested on *nix based systems.

==
?- gbn(debug(true)).
% Turning debugging on for predicate handle: gbn(gbn)
% Options: [$restore(gbn,debug,false),copy(false),data(pack(gbn/data/asia.dat)),display_dot(svg),odir(_33202),std_output(std_file)]
% Output directory: 'asia-20.09.27'
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

@author nicos angelopoulos
@see gbn_version/2

*/

/** gbn_module.

Documentation predicate to give anchor to the overall module documentation.

*/
gbn_module.
