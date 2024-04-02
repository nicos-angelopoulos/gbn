
gbn_mtx_subs_defaults( [res(Res)] ) :-
    os_unique( res, Res, create(false) ).

/** gbn_mtx_subs( +Mtx, +Pred, +Opts ).

Run gbn on a number of subsets of Mtx as defined by Pred/+2
Pred should be defined with additional 2 arguments, first for the name/identifier
of the subset and second for list of pairs of the form: Column-Condition
which be passed as 2nd and 3rd args of mtx_subset/4.

*/
gbn_mtx_subs( MtxIn, Pred, Args ) :-
    options_append( gbn_mtx_subs, Args, Opts ),
    options( to_bn_mtx(BnMtxG), Opts ),
    options( res(Res), Opts ),
    os_make_path( Res ),
    mtx( MtxIn, Mtx ),
    call( Pred, Name, Conds ),
    gbn_mtx_conditions( Conds, Mtx, Contx ),
    call( BnMtxG, Contx, BnMtx ),
    os_dir_stem_ext( data, Name, csv, CsvF ),
    mtx( CsvF, BnMtx ),
    gbn_mtx_dat( CsvF, DatF ),
    debuc( gbn(mtx_subs), 'DatF: ~p', DatF ),
    numlist( 1, 10, ToTen ),
    os_path( Res, Name, ThisRes ),
    gbn( [dir(ThisRes),data(DatF),multiple(edge_penalty,ToTen)] ),
    fail.
gbn_mtx_subs( _MtxIn, _Pred, _Args ).

gbn_mtx_conditions( [], Mtx, Mtx ).
gbn_mtx_conditions( [Cid-Goal|T], Mtx, Contx ) :-
    mtx_subset( Mtx, Cid, Goal, Mid ),
    gbn_mtx_conditions( T, Mid, Contx ).
