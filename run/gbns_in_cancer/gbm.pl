
:- use_module(library(lib)).
:- use_module(library(debug)).

:- lib(gbn).
:- lib(debug_call).

:- debug(gbm).

gbm :-
    debug_call( gbm, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/gbm/gbm_min05.dat'),
    E = 1,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debug_call( gbm, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debug_call( gbm, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debug_call( gbm, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debug_call( gbm, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debug_call( gbm, end, true ).
    
