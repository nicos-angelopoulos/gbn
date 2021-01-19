
:- lib(mtx).

/** gbn_mtx_paired( +Mtx, +Cid1, +Cid2, -PAMtx ).

Create a presence/absence matrix from two colums in Mtx. 

For each value of Cid1 a row is created in PAMtx.
For each value in Cid2 a column is created in PAMtx.

Values in PAMtx are binary and show presence (1) or absence of the
particular Cid1/Cid2 values in Mtx.

PAMtx contains a header (all the Cid2 values, ordered) but now row names.

==
?- Mtx = [r(a,b,c),r(1,2,1),r(2,2,2)],
   gbn_mtx_paired( Mtx, a, c, PA ).

PA = [row(1, 2), row(1, 0), row(0, 1)].

?- Mtx = [r(ps,xs,ys),r(p1,x1,y1),r(p2,x2,y2),r(p1,x2,y3),r(p3,x3,y3)],
   gbn_mtx_paired( Mtx, ps, xs, PA ).

PA = [row(x1, x2, x3), row(1, 1, 0), row(0, 1, 0), row(0, 0, 1)].

?- Mtx = [r(ps,xs,ys),r(p1,x1,y1),r(p2,x2,y2),r(p1,x2,y3),r(p3,x3,y3),r(p1,x1,y1)],
   gbn_mtx_paired( Mtx, ps, xs, PA ).

Mtx = [r(ps, xs, ys), r(p1, x1, y1), r(p2, x2, y2), r(p1, x2, y3), r(p3, x3, y3), r(p1, x1, y1)],
PA = [row(x1, x2, x3), row(1, 1, 0), row(0, 1, 0), row(0, 0, 1)].

==
*/
gbn_mtx_paired( MtxIn, Cid1, Cid2, PAM ) :-
    mtx( MtxIn, Mtx ),
    mtx_column_set( Mtx, Cid1, Set1 ),
    mtx_column_set( Mtx, Cid2, Set2 ),
    mtx_columns_kv( Mtx, Cid1, Cid2, KVs, _Cnms, _Cpos),
    sort( KVs, KVo ),
    gbn_columns_pairs_pam( KVo, Set1, Set2, PAMRows ),
    PAMHdr =.. [row|Set2],
    PAM = [PAMHdr|PAMRows]. 

gbn_columns_pairs_pam( [], Ks, _, [] ) :-
    ( Ks == [] -> true; trow( residual_ks_in_gbn_mtx_paired(Ks) ) ).
gbn_columns_pairs_pam( [K-V|KVs], [K|Ks], Vset, [KRow|Rows] ) :-
    !,
    gbn_set_value_pa_row( Vset, V, RVset, PAList, PATail ),
    gbn_columns_pairs_pam_k( KVs, K, RVset, RemKVs, PATail ),
    KRow =.. [row|PAList],
    gbn_columns_pairs_pam( RemKVs, Ks, Vset, Rows ).
gbn_columns_pairs_pam( [K-V|_KVs], [K1|Ks], _Vset, _ ) :-
    throw( gbn_mtx_paired_key_mismatch(K,K1,V,Ks) ).

gbn_set_value_pa_row( Vset, V, RemVset, PAList, PATail ) :-
    gbn_set_value_pa_row_1( Vset, V, RemVset, PAList, PATail ),
    !.
gbn_set_value_pa_row( [], V, _RemVset, _PAList, _PATail ) :-
    throw( gbn_mtx_paired_unable_to_establish_presence_of(V) ).
    
gbn_set_value_pa_row_1( [V|Vs], V, RVset, PAList, PATail ) :-
    !,
    PAList = [1|PATail],
    RVset = Vs.
gbn_set_value_pa_row_1( [_V1|Vs], V, RVset, PAList, PATail ) :-
    PAList = [0|PACont],
    gbn_set_value_pa_row_1( Vs, V, RVset, PACont, PATail ).

gbn_columns_pairs_pam_k( [], _K, Vset, [], PAList ) :-
    findall( 0, member(_,Vset), PAList ).
gbn_columns_pairs_pam_k( [K-V|T], K, Vset, RemKVs, PAList ) :-
    !,
    gbn_set_value_pa_row( Vset, V, RVset, PAList, PACont ),
    gbn_columns_pairs_pam_k( T, K, RVset, RemKVs, PACont ).
gbn_columns_pairs_pam_k( [H|T], _K, Vset, RemKVs, PAList ) :-
    RemKVs = [H|T],
    findall( 0, member(_,Vset), PAList ).
