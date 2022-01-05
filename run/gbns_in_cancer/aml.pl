
:- use_module(library(lib)).

:- lib(gbn).
:- lib(debug_call).

:- debuc(aml).

/** aml.

Run GBN on the AML cancer dataset with \epsilon = 7 producing all outputs
as shown on GBN in cancer paper.

==
?- [cancer(aml)].
?- aml.
% Starting: aml
% Starting: fisher_nets
% Starting: fam_hmaps
% Starting: gates_nets
% Starting: svg_legend
% Finished: aml
true.
==

*/

aml :-
    debuc( aml, start, true ),
    % DatF = 'data/aml/aml_min60.dat',
    DatF = pack('gbn/data/gbns_in_cancer/aml/aml_min60.dat'),
    E = 7,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debuc( aml, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debuc( aml, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debuc( aml, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debuc( aml, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debuc( aml, end, true ).
    
/** aml_full.

Run GBN on the aml cancer dataset for \epsilon =1...10.

*/
aml_full :-
    debuc( aml, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/aml/aml_min60.dat'),
    numlist( 1, 10, Es ),
    GBNOpts = [ data(DatF),
                multiple(edge_penalty,Es),
                dir(OsOdir)
              ],
    gbn( GBNOpts ),
        debuc( aml, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debuc( aml, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debuc( aml, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debuc( aml, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debuc( aml, end, true ).
    
