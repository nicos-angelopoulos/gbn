
:- use_module(library(lib)).
:- use_module(library(debug)).

:- lib(gbn).
:- lib(debug_call).

:- debug(mpn).

mpn :-
    debug_call( mpn, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/mpn/mpn_min05.dat'),
    E = 3,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debug_call( mpn, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debug_call( mpn, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debug_call( mpn, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debug_call( mpn, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debug_call( mpn, end, true ).
    
