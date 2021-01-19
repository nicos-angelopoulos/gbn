
:- use_module(library(lib)).
:- use_module(library(debug)).

:- lib(gbn).
:- lib(debug_call).

:- debug(coa).

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
    
