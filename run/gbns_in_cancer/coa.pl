
:- use_module(library(lib)).
:- use_module(library(debug)).

:- lib(gbn).
:- lib(debug_call).

:- debug(coa).

/** coa

Run GBN on the colorectal cancer dataset with \epsilon = 1 producing all outputs
as shown on GBN in cancer paper.

*/

coa :-
    debug_call( coa, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/coa/coa_min05.dat'),
    E = 1,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debug_call( coa, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debug_call( coa, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debug_call( coa, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debug_call( coa, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debug_call( coa, end, true ).
    
/** coa_all.

Run GBN on the colorectal cancer dataset for \epsilon =1...10.

*/
coa_all :-
    debug_call( coa, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/coa/coa_min05.dat'),
    numlist( 1, 10, Es ),
    GBNOpts = [ data(DatF),
                multiple(edge_penalty,Es),
                dir(OsOdir)
              ],
    gbn( GBNOpts ),
        debug_call( coa, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debug_call( coa, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debug_call( coa, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debug_call( coa, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debug_call( coa, end, true ).
    
