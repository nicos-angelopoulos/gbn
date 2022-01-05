% :- ensure_loaded(aml_plots).
:- use_module(library(apply)).
:- use_module(library(lists)).
:- use_module(library(lib)).

:- lib(mtx).
:- lib(real).
:- lib(b_real).
:- lib(os_lib).
:- lib(stoics_lib).
:- lib(debug_call).

:- debuc(pplots).

pop_plots :-
    debuc( pplots, start, true ),
    @ mkdir(tmp_pop_plots),  % will fail and execution stopped if this dir exists
    maplist( pop_dset, [coa,gbm,mpn,mye] ),
    debuc( pplots, end, true ).

pop_dset( Dset ) :-
    debuc( pplots, start, Dset ),
    os_ext( dat, Dset, Dfile ),
    at_con( ['../../../data/gbns_in_cancer',Dset,Dfile], '/', Path ),
    mtx( Path, Mtx, sep(0' ) ),
    debuc( pplots, dims, mpn/Mtx ),
    mtx_lists( Mtx, Lists ),
    maplist( min_ones(Lists), [0,20,40,60,80], Pops ),
    debuc( pplots, 'KVs: ~w', [Pops] ),
    FClrs = ["#779ECB"],
    pop_name( Dset, SetName ),
    at_con( ['Effect of cut-off on number of variables (',SetName,').'], '', Title ),
    Ppts = [ flip(false), 
             geom_bar(empty),
             labels('Cut-off for minimum occurances','Number of Variables',Title),
             fill_colours(FClrs),
             panel_theme(axes)
           ],
    gg_bar_plot( Pops, Ppts ),
    os_ext( pdf, Dset, PdfSet ),
    os_postfix( [min,pop], PdfSet, PdfF ),
    os_path( tmp_pop_plots, PdfF, PdfP ),
    <- ggsave( +PdfP ),
    debuc( pplots, end, Dset ).

min_ones( MtxLists, AtLeast, AtLeast-Clms ) :-
    include( at_least_ones(AtLeast), MtxLists, RemLists ),
    length( RemLists, Clms ).

at_least_ones( AtLeast, [_,_|Nums] ) :-
    sum_list( Nums, Sum ),
    AtLeast =< Sum.

pop_name( mye, myeloma ).
pop_name( mpn, 'MPN' ).
pop_name( coa, 'colorectal' ).
pop_name( gbm, 'glioblastoma' ).
