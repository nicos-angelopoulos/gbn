
:- use_module(library(lib)).
:- use_module(library(debug)).

:- lib(gbn).
:- lib(debug_call).

:- debug(aml).

aml :-
    debug_call( aml, start, true ),
    DatF = 'data/aml/aml_min60.dat',
    E = 7,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debug_call( aml, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debug_call( aml, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debug_call( aml, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debug_call( aml, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debug_call( aml, end, true ).
    
