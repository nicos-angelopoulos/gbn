
% :- lib(suggests(svg)).
% :- lib(suggests(mtx)).
% :- lib(suggests(real)).

% :- lib(stoics_lib:kv_decompose/3).
% :- lib(stoics_lib:list_proportions/3).

gbn_svg_legend_defaults( Defs ) :-
    Defs = [
                adds_lines(gbn:gbn_svg_add_lines),
                theme(stack),
                placement(bottom),
                has_postfix(gated),
                has_postfix(fclr),
                has_postfix(fisher)
        ].

/**  gbn_svg_legend(Opts)

A simple gbn shell to svg_legend/1.

Opts
  See svg_legend/1.

==
==

@author nicos angelopoulos
@version  0.1 2018/2/26
@see svg_legend/1

*/
gbn_svg_legend( Args ) :-
    options_append( gbn_svg_legend, Args, Opts ),
    svg_legend( Opts ).

gbn_svg_add_lines( stack, Xi, Yi, File, Add, MaX, MaY ) :-
    ( atom_concat(Stem,'fclr_gated.svg',File) ->
        atom_concat(Stem,'gates_best.csv',GatesF),
        mtx( GatesF, GatesMtx, [convert(true),match(false)] ),
        findall( Gate, (  member(Row,GatesMtx), arg(6,Row,RowGatesAtm),
                          atom_to_term(RowGatesAtm,RowGates,[]),
                          gates_term_functor(RowGates,Gate)
                       )
                        , GatesNestBag ),
        flatten( GatesNestBag, GatesBag ),
        gates_sort( GatesBag, Gates ),
        ( Gates == [] ->
            Y5 = Yi, GatesAtms = []
            ;
            svg_add_lines_gates( Gates, Xi, Yi, Y5, GatesAtms )
        )
        ;
        Y5 = Yi, GatesAtms = []
    ),
    Xt is Xi + 30,
    svg_add_lines_compos( [node,pval,cooc,muex,slide], File, Xi, Xt, Y5, 0, -, MaX/MaY, CompoLnsNest ),
    flatten( CompoLnsNest, CompoLns ),
    append( GatesAtms, CompoLns, AllLns ),
    % maplist( atom_codes, GatesAtms, Add ).
    maplist( atom_codes, AllLns, Add ).

svg_add_lines_gates( [], _Xi, Y, Y, [] ).
svg_add_lines_gates( [G|Gs], Xi, Yi, Y5, [Ln1,Ln2,Ln3,Ln4|GatesAtms] ) :-
    svg_add_lines_gate( G, Xi, Yi, Ln1, Ln2, Ln3, Ln4 ),
    Yj is Yi - 20,
    svg_add_lines_gates( Gs, Xi, Yj, Y5, GatesAtms ).

% svg_add_lines_gate( xor, _Xi, _Yi, _Yj, _L1, _L2, _L3, _L4 ) :-
    % throw( xor(ing) ).
svg_add_lines_gate( xor, Xi, Yi, Lat, Mat, Nat, Oat ) :-
    % OR GATE
    T1 = '<text text-anchor="start" x="',
    T2 = '" y="',
    T3 = '" font-family="Times,serif" font-size="14.00">', 
    T4 = '</text>',
    OrX is Xi + 30,
    OrY is Yi,
    OrTxt = 'XOR gate',
    Mat = ' <g id="legXOR" class="node"><title>LegXOR</title>',
    N1 = ' <polygon fill="#87a96b" stroke="#87a96b" points="',
    maplist( plus(Xi), [13,15,19,21,13], Xs ),
    % maplist( plus(Yi), [-54,-51,-54,-58,-58,-54], Ys ),
    % maplist( plus(Yi), [74,78,74,71,71,74], Ys ),
    maplist( plus(Yi), [-8,-1,-1,-8,-8], Ys ),
    % Whatevs = '6,-23 10,-20 14,-23 14,-27 6,-27 6,-23',
    maplist( atom_concat(','), Ys, ComYs ),
    maplist( atom_concat, Xs, ComYs, XYs ),
    atomic_list_concat( XYs, ' ', Poly ),
    N2 = '"/>',
    Oat  = '</g>',
    atomic_list_concat( [T1,OrX,T2,OrY,T3,OrTxt,T4], Lat ),
    atomic_list_concat( [N1,Poly,N2], Nat ).

svg_add_lines_gate( not, Xi, Yi, NotL1, NotL2, NotL3, NotL4 ) :-
    % NOT gate 
    NotTxX is Xi + 30,
    NotTxY is Yi,
    NotX is Xi + 17,
    NotY is Yi - 5,
    NotL1a = '<text text-anchor="start" x="',
    NotL1b = '" y="',
    NotL1c = '" font-family="Times,serif" font-size="14.00">NOT gate</text>',
    NotL2 = '<g id="legNot" class="node"><title>LegNot</title>',
    NotL3a = '<ellipse fill="#87a96b" stroke="#87a96b" cx="', % asparagus
    NotL3b = '" cy="',
    NotL3c = '" rx="3" ry="3"/>',
    NotL4 = '</g>',
    atomic_list_concat( [NotL1a,NotTxX,NotL1b,NotTxY,NotL1c], NotL1 ),
    atomic_list_concat( [NotL3a,NotX,NotL3b,NotY,NotL3c], NotL3 ).
svg_add_lines_gate( and, Xi, Yi, AndLat, AndMat, AndNat, Aat ) :-
    % AND GATE
    T1 = '<text text-anchor="start" x="',
    T2 = '" y="',
    T3 = '" font-family="Times,serif" font-size="14.00">', 
    T4 = '</text>',
    AndX is Xi + 30,
    AndY is Yi,
    AndTxt = 'AND gate',
    AndMat = ' <g id="legAND" class="node"><title>LegAND</title>',
    AndN1 = ' <polygon fill="#87a96b" stroke="#87a96b" points="',
    maplist( plus(Xi), [13,17,21,13], AndXs ),
    maplist( plus(Yi), [-8,-1,-8,-8], AndYs ),
    % Whatevs = '6,-23 10,-20 14,-23 14,-27 6,-27 6,-23',
    maplist( atom_concat(','), AndYs, AndComYs ),
    maplist( atom_concat, AndXs, AndComYs, AndXYs ),
    atomic_list_concat( AndXYs, ' ', AndPoly ),
    AndN2 = '"/>',
    Aat  = '</g>',
    atomic_list_concat( [T1,AndX,T2,AndY,T3,AndTxt,T4], AndLat ),
    atomic_list_concat( [AndN1,AndPoly,AndN2], AndNat ).
svg_add_lines_gate( or, Xi, Yi, Lat, Mat, Nat, Oat ) :-
    % OR GATE
    T1 = '<text text-anchor="start" x="',
    T2 = '" y="',
    T3 = '" font-family="Times,serif" font-size="14.00">', 
    T4 = '</text>',
    OrX is Xi + 30,
    OrY is Yi,
    OrTxt = 'OR gate',
    Mat = ' <g id="legOR" class="node"><title>LegOR</title>',
    N1 = ' <polygon fill="#87a96b" stroke="#87a96b" points="',
    maplist( plus(Xi), [13,17,21,21,13,13], Xs ),
    % maplist( plus(Yi), [-54,-51,-54,-58,-58,-54], Ys ),
    % maplist( plus(Yi), [74,78,74,71,71,74], Ys ),
    maplist( plus(Yi), [-5,-1,-5,-8,-8,-5], Ys ),
    % Whatevs = '6,-23 10,-20 14,-23 14,-27 6,-27 6,-23',
    maplist( atom_concat(','), Ys, ComYs ),
    maplist( atom_concat, Xs, ComYs, XYs ),
    atomic_list_concat( XYs, ' ', Poly ),
    N2 = '"/>',
    Oat  = '</g>',
    atomic_list_concat( [T1,OrX,T2,OrY,T3,OrTxt,T4], Lat ),
    atomic_list_concat( [N1,Poly,N2], Nat ).

gates_sort( GatesBag, Gates ) :-
    sort( GatesBag, GatesOrd ),
    reverse( GatesOrd, GatesInRev ),
    ( select(not,GatesInRev,GatesNoNot) ->
        Gates = [not|GatesNoNot]
        ;
        Gates = GatesInRev
    ).

gates_term_functor( Atom, Funcs ) :-
    atomic( Atom ),
    !,
    Funcs = [].
gates_term_functor( Term, [Token|FuncNest] ) :-
    Term =.. [Func|Args],
    gates_term_functor_atom( Func, Token ),
    maplist( gates_term_functor, Args, FuncNest ).

gates_term_functor_atom( o, or  ) :- !.
gates_term_functor_atom( a, and ) :- !.
gates_term_functor_atom( n, not ) :- !.
gates_term_functor_atom( x, xor ) :- !.
    % throw( xor(ing) ).
gates_term_functor_atom( Oth, _ ) :- !,
    throw( unknown_gate_functor(Oth) ).

gbn_svg_line_segs_as( A1, A2, A3, A4 ) :-
    A1 = '<text text-anchor="start" x="', % 856
    A2 = '" y="', % -85
    A3 = '" font-family="Times,serif" font-size="14.00">',
    A4 = '</text>'.
gbn_svg_line_segs_bs( B1, B2, B3 ) :-
    B1 = '<path fill="none" stroke-width="2.8" stroke="#',
    B2 = '" d="M',
    B3 = ' l 25,0"/>'.
 
svg_add_lines_compos( [], _File, _Xi, _Xt, Yi, MaX, _Dir, MaX/Yi, [] ).
svg_add_lines_compos( [C|Cs], File, Xi, Xt, Yi, CurMaX, Dir, Max, [Lns|TLns] ) :-
    debuc( gbn(gbn), 'Adding compo: ~w', C ),
    ( svg_add_lines_compo( C, File, Xi, Xt, Yi, CurMaX, IncY, NxtMaX, Lns ) ->
        true
        ;
        debuc( gbn(gbn), 'failed compo: ~w', C ),
        Lns = [],
        IncY is 0,
        NxtMaX is CurMaX
    ),
    Expr =.. [Dir,Yi,IncY],
    Yj is Expr,
    % Yj is Yi + IncY,
    svg_add_lines_compos( Cs, File, Xi, Xt, Yj, NxtMaX, Dir, Max, TLns ).

svg_add_lines_compo( slide, _File, Xi, Xt, Yi, CurMaX, 20, NxtMaX, [I,J] ) :-
    % compo: slide
    gbn_svg_line_segs_as( A1, A2, A3, A4 ),
    L5 = 'Fisher test odds',
    atomic_list_concat( [A1,Xt,A2,Yi,A3,L5,A4], '', I ),
    J1 = '<polygon points="',
    Tx1 is Xi, Ty1 is Yi,
    Tx2 is Tx1 + 25, Ty2 is Ty1,
    Tx3 is Tx2, Ty3 is Ty2 - 6,
    J1b = '" style="fill:#BF812D;stroke:#35978F;stroke-width:1.4" />',
    NxtMaX is max(CurMaX,200),  % fixme: 
    atomic_list_concat( [J1,Tx1,',',Ty1,' ',Tx2,',',Ty2,' ',Tx3,',',Ty3,J1b], J ).

svg_add_lines_compo( cooc, _File, Xi, Xt, Yi, CurMaX, 20, NxtMaX, [A,E] ) :-
    gbn_svg_line_segs_as( A1, A2, A3, A4 ),
    gbn_svg_line_segs_bs( B1, B2, B3 ),
    L1 = 'Co-occur (shown odds=4)',
    atomic_list_concat( [A1,Xt,A2,Yi,A3,L1,A4], '', A ),
    Ya is Yi - 5,
    NxtMaX is max(CurMaX,200),  % fixme:
    atomic_list_concat( [B1,'35978F',B2,Xi,',',Ya,B3], E ).

svg_add_lines_compo( muex, _File, Xi, Xt, Yi, CurMaX, 20, NxtMaX, [A,E] ) :-
    gbn_svg_line_segs_as( A1, A2, A3, A4 ),
    gbn_svg_line_segs_bs( B1, B2, B3 ),
    L2 = 'Mut.excl (shown odds=0.25)',
    Ya is Yi - 5,
    atomic_list_concat( [A1,Xt,A2,Yi,A3,L2,A4], '', A ),
    NxtMaX is max(CurMaX,200),  % fixme:
    atomic_list_concat( [B1,'BF812D',B2,Xi,',',Ya,B3], E ).

svg_add_lines_compo( pval, File, Xi, Xt, Yi, CurMaX, Pad, NxtMaX, Lns ) :-
    % ( File == 'dutch_driver_muts_min20_fisher.svg' -> trace; true ),
    gbn_svg_line_segs_as( A1, A2, A3, A4 ),
    gbn_svg_line_segs_bs( _B1, B2, _B3 ),
    os_ext( _, csv, File, CsvFPrv ),
    ( exists_file(CsvFPrv) -> CsvF = CsvFPrv
            ; 
            ( os_postfix(gated,CsvF,CsvFPrv) -> true; CsvF=CsvFPrv)
    ),
    ( (exists_file(CsvF),mtx(CsvF,Mets,convert(true))) ->
        maplist( arg(3), Mets, Pvals ),
        max_list( Pvals, MaxPv )
        ; 
        % if there is not 
        MaxPv is 0
    ),
    % ( mtx(CsvF,Mets) -> true; throw( could_not_mtx_csv_file(CsvF) ) ),
    ( MaxPv < 0.05  ->   % fixme: = means no significance...
        Pad is 0,
        Lns = [],
        NxtMaX is CurMaX
        ;
        NxtMaX is max(CurMaX,200),  % fixme:
        Pad is 20,
        L3 = '0.05 &lt; q.val',
        D1 = '<path fill="none" stroke-width="2.8" stroke-dasharray="5,2" stroke="#',
        atomic_list_concat( [A1,Xt,A2,Yi,A3,L3,A4], '', C ),
        Yim is Yi - 5,
        C3 = ' l 11,0"/>',
        atomic_list_concat( [D1,'35978F',B2,Xi,',',Yim,C3], DA ),
        Xib is Xi + 14,
        atomic_list_concat( [D1,'BF812D',B2,Xib,',',Yim,C3], DB ),
        Lns = [C,DA,DB]
    ).

svg_add_lines_compo( node, _File, Xi, Xt, Yi, CurMaX, 20, NxtMaX, [G,H] ) :-
    gbn_svg_line_segs_as( A1, A2, A3, A4 ),
    gbn_svg_line_segs_bs( _B1, B2, B3 ),
    L4 = '# events (shown med=', L3b = ',max=', L3c = ')',
    % svg_dat_file( File, DatFile ),
    gbn_res_dir_dat_file( '.', DatFile ),
    debuc( gbn(gbn), 'dat.file(~p).', DatFile ),
    mtx( DatFile, DatMtx, [sep(0' ),convert(true)] ),
    mtx_value_column_frequencies( DatMtx, 1, Freqs ),
    kv_decompose( Freqs, _Lbls, Times ),
    TmMedian <- as.integer( median( Times ) ),
    TmMax <- as.integer( max( Times ) ),
    list_proportions( Times, Propos, to_range(r(1,4)) ),
    NwMedian <- median( Propos ),
    % Yl is Yk + 20,
    atomic_list_concat( [A1,Xt,A2,Yi,A3,L4,TmMedian,L3b,TmMax,L3c,A4], '', G ),
    G1 = '<path fill="none" stroke-width="', G1b= '" stroke="#',
    Yg is Yi - 5,
    NxtMaX is max(CurMaX,230),  % fixme:
    atomic_list_concat( [G1,NwMedian,G1b,'000000',B2,Xi,',',Yg,B3], H ).

/*
svg_dat_file( File, DatFile ) :-
    atomic_list_concat( [FStem|_], '-', File ),
    os_ext( dat, FStem, DatFile ),
    exists_file( DatFile ),
    !.
svg_dat_file( File, DatFile ) :-
    os_postfix( fisher, Base, File ),
    os_ext( _OldExt, dat, Base, DatFile ),
    exists_file( DatFile ),
    !.
svg_dat_file( File, DatFile ) :-
    atomic_list_concat( [FStem|_], '_', File ),
    os_ext( dat, FStem, DatFile ),
    exists_file( DatFile ),
    !.
svg_dat_file( File, _DatFile ) :-
    throw( cannot_locate_dat_file_for(File) ).
    */
