
% :- lib(real).
% :- lib(stoics_lib:skim/3).
% :- lib(stoics_lib:at_con/3).

gbn_family_gates_defaults( [simplify(aox),xor(true),not(false)] ).


/** gbn_family_gates( +Child, +Parents, +Mtx, -Gatrix, +Opts ).

    For a Child and list of Parents variables that are columns in Mtx, <br>
    Gatrix is the matrix formed of the columns described below.        <br>
    Note that values in the Child and Parents columns should be binary for the predicate to work.<br>

Opts 
  * not(Not=false)
    whether to include not gates
  * simplify(Sify)
    defines the type of simplication on the gates term (_aox_ default, else no simplification)
  * xor(Xor=true)
    whether to include xor gates

Gatrix
  * test_statistic
     the result of Fisher test on the logic formula
  * complexity
     integer- number of logic gates in formula
  * direction
     mut_excl, or co_occur
  * odds
     value on which the direction of relation between the formula and Child was derived
  * pval
     pvalue
  * formula
     logical gates formula

The first row holds the Child and the Parents consecutively.  <br>
Gatrix is passed through mtx/2 before being returned. Mtx is similarly treated as mtx/2 input.

@author nicos angelopoulos
@version  0.2 added options and ability to exclude not (currently false is the default)

*/
gbn_family_gates( Child, Parents, MtxIn, Gatrix, Args ) :-
    options_append( gbn_family_gates, Args, Opts ),
    options( [simplify(Sify),xor(Xor),not(Not)], Opts ),
    gbn_family_gates( Child, Parents, MtxIn, Xor, Not, Sify, Gatrix ).

gbn_family_gates( Child, Parents, MtxIn, Wxor, Not, Simplify, Gatrix ) :-
    gates( Parents, Wxor, Not, Simplify, Gates ),
    mtx_columns_values( MtxIn, Mtx, header_pair(true) ),
    memberchk( Child-ChildValues, Mtx),
    findall( Pa-PaValues, (member(Pa,Parents),memberchk(Pa-PaValues,Mtx)), VVs ), 
    childvs <- ChildValues,
    gate_expr_pvals( Gates, VVs, childvs, Pates ),
    sort( Pates, Patrix ),
    gbn_expr_pvals_epsilon_correct( Patrix, Catrix ),
    sort( Catrix, Gatrix ).

gbn_expr_pvals_epsilon_correct( [row(P,C,D,E,A,M)|T], Corrected ) :-
    gbn_expr_pvals_epsilon_correct_1( [row(P,C,D,P,E,A,M)|T], Corrected ).
    

gbn_expr_pvals_epsilon_correct_1( [Row], [Row] ) :- !.
gbn_expr_pvals_epsilon_correct_1( [row(P1,C1,D1,N1,E1,A1,M1),row(P2,C2,D2,E2,A2,M2)|T], [CRow|Catrix] ) :-
    !,
    CRow = row(P1,C1,D1,N1,E1,A1,M1),
    ( catch(log10(P1) + 0.00001 < log10(P2),_,fail) ->  N2 = P2; N2 = P1 ),
    % ( P1 + epsilon < P2 -> N2 = P2; N2 = P1 ),
    gbn_expr_pvals_epsilon_correct_1( [row(N2,C2,D2,P2,E2,A2,M2)|T], Catrix ).

gate_expr_pvals( [], _VVs, _ChildValues, [] ).
gate_expr_pvals( [G|Gs], VVs, ChildValues, [Row|Prs] ) :-
    % debug( gbn(fisher_test), 'gate: ~w', [G] ),
    % debug( gbn(fisher_test), 'vvs: ~w', [VVs] ),
    G =.. [Gn|Gargs],
    % gate_eval( G, VVs, Cx, Vect ),
    gate_eval( Gn, Gargs, VVs, ChildValues, Cx, Vect, SubMts ),
    gate_vector_pairwise_metrics( Vect, ChildValues, Pval, Odds, Dir ),
    % G = _-Formula,
    term_to_atom( G, FormAtm ),
    % Row = row(Pval,Cx,Dir,Odds,FormAtm,SubMts),
    % here() too:
    Row = row(Pval,Cx,Dir,Odds,FormAtm,SubMts),
    gate_expr_pvals( Gs, VVs, ChildValues, Prs ).

gate_vector_pairwise_metrics( Vect, ChildValues, Pv, Pe, Dir ) :-
    ( (\+ sort(Vect,[_]),catch(ft <- fisher.test( Vect, ChildValues),_,fail)) ->
        % debug( gbn(fisher_test), 'Run fisher.test on vect: ~w', [Vect] ),
        % ChildValuesL <- ChildValues,
        % debug( gbn(fisher_test), 'Run fisher.test against: ~w', [ChildValuesL] ),
        PvPrv <- ft$p.value,
        PePrv <- ft$estimate,
        ( number(PePrv) -> 
            Pv is PvPrv, Pe is PePrv,
            ( Pe > 1 -> Dir = co_occur; Dir = mut_excl )
            ;
            Pv = 1, Pe = 1,
            Dir = null
        )
        ;
        Pv = 1,
        Pe = 1,
        Dir = null
    ).

/** gates( +Vars, +WxorB, +NotB, +Simplify, +GateExprs ).
    
Create complete set of gate expression that involve all of elements in Vars.<br>
WxorB is the boolean controlling inclusion of Xor gate.<br>
Simplify is the token controlling simplification regime,<br>
_aox_ simplification collapses all binaries and any other token performs no simplification.

without xor: 2:16, 3:96, 4:1024, 5:8192, 6:65536, 7:524288, out of global<br>
with xor: 2:24, 3:288, 4:3456, 5:41472, 6:497664, 7:out of global...

==
?- gbn:gates( [b,c], true, false, Gxs ), maplist( writeln, Gxs ), length( Gxs, Len ).
a(b,c)
n(a(b,c))
o(b,c)
n(o(b,c))
...
x(n(b),n(c))
n(x(n(b),n(c)))
Gxs = [a(b, c), n(a(b, c)), o(b, c), n(o(b, c)), x(b, c), n(x(b, c)), a(b, n(c)), n(a(..., ...)), o(..., ...)|...],
Len = 24.

?- gbn:gates( [b,c], false, false, Gxs ), maplist( writeln, Gxs ), length( Gxs, Len ).
a(b,c)
n(a(b,c))
o(b,c)
n(o(b,c))
a(b,n(c))
n(a(b,n(c)))
o(b,n(c))
n(o(b,n(c)))
a(n(b),c)
n(a(n(b),c))
o(n(b),c)
n(o(n(b),c))
a(n(b),n(c))
n(a(n(b),n(c)))
o(n(b),n(c))
n(o(n(b),n(c)))
Gxs = [a(b, c), n(a(b, c)), o(b, c), n(o(b, c)), a(b, n(c)), n(a(b, n(c))), o(b, n(c)), n(o(..., ...)), a(..., ...)|...],
Len = 16.

?- gbn:gates( [b,c,d], true, false, Gxs ), maplist( writeln, Gxs ), length( Gxs, Len ).
a(a(b,c),d)
n(a(a(b,c),d))
o(a(b,c),d)
...
Len = 288.
==

?- gbn:gates( [b,c,d], true, aox, Gxs ), maplist( writeln, Gxs ), length( Gxs, Len ).
a(b,c,d)
n(a(b,c,d))
o(a(b,c),d)
n(o(a(b,c),d))
....
Len = 288.

*/
gates( [], _Wxor, _Not, _Simfy, [] ).
gates( [A,B|T], Wxor, Not, Simfy, Gates ) :-
    gates( [A], Wxor, Not, Simfy, Ags ),
    gates( [B], Wxor, Not, Simfy, Bgs ),
    findall( ABinGs, (member(AnAg,Ags),member(ABg,Bgs),gates_binary(AnAg,ABg,Wxor,Simfy,ABinGs)), BinGsNest ),
    flatten( BinGsNest, BinGs ),
    !,
    findall( ATGs, (member(BinG,BinGs),gates([BinG|T],Wxor,Not,Simfy,ATGs)), TGs ),
    flatten( TGs, Gates ).
gates( [A], _Wxor, Not, _Simfy, Gates ) :-
    gates_not( Not, A, Gates ).

gates_not( true, A, [A,n(A)] ).
gates_not( false, A, [A] ).

% gates_binary( A, B, [a(A,B),o(A,B)] ).
% gates_binary( A, B, [a(A,B),o(A,B),x(A,B)] ).
gates_binary( A, B, Wxor, Simfy, Gates ) :-
    gates_binary_with_xor( Wxor, A, B, Binaries ), % fixme: protect
    maplist( gate_binary_simplify(Simfy), Binaries, Gates ).

gates_binary_with_xor( true, A, B, [a(A,B),o(A,B),x(A,B)] ).
gates_binary_with_xor( false, A, B, [a(A,B),o(A,B)] ).

gate_binary_simplify( aox, Gate, Simfied ) :- !,
    functor( Gate, Gn, 2 ),
    gate_binary_arg_gate( Gate, 1, Arg1, A1n, A1Args ),
    gate_binary_arg_gate( Gate, 2, Arg2, A2n, A2Args ),
    gate_binary_simply_aox( Gn, A1n, A2n, Arg1, Arg2, A1Args, A2Args, Simfied ).
gate_binary_simplify( _, Gate, Gate ).   % makes it the default to any other value
    
gate_binary_arg_gate( Gate, Pos, Arg, An, AArgs ) :-
    arg( Pos, Gate, Arg ),
    Arg =.. [An|AArgs].

gate_binary_simply_aox( Gn, Gn, Gn, _Arg1, _Arg2, A1Args, A2Args, Simfied ) :- 
    !,
    append( A1Args, A2Args, Args ),
    Simfied =.. [Gn|Args].
gate_binary_simply_aox( Gn, Gn, _A2n, _Arg1, Arg2, A1Args, _A2Args, Simfied ) :- 
    !,
    append( A1Args, [Arg2], Args ),
    Simfied =.. [Gn|Args].
gate_binary_simply_aox( Gn, _A1n, Gn, Arg1, _Arg2, _A1Args, A2Args, Simfied ) :- 
    !,
    Simfied =.. [Gn,Arg1|A2Args].
gate_binary_simply_aox( Gn, _A1n, _A2n, Arg1, Arg2, _A1Args, _A2Args, Simfied ) :- 
    Simfied =.. [Gn,Arg1,Arg2].

/* 18.02.13: i think the following is incorrect:
gates( Parts, Gates ) :-
    findall( [P]-P, member(P,Parts), PairParts ),
    gates( PairParts, [], [], Gates ).

gates( [], Seen, Acc, Gates ) :-
    gates_cont( Acc, Seen, Gates ).
gates( [P|Parts], Seen, Acc, Gates ) :-
    ( ord_memberchk(P,Seen) -> 
        Qarts = Parts,
        Neen = Seen,
        Nxt = Acc,
        Gates = Tates
        ;
        gate_expand( Parts, P, Parted, Tarted ),
        gate_expand( Seen, P, Tarted, [] ),
        gates_acc( Parted, P, Acc, Nxt, Gates, Tates ),
        P = Set-Expr,
        ( Expr = n(_) -> Parts = Qarts; Qarts = [Set-n(Expr)|Parts] ),
        ord_add_element( Seen, P, Neen )
    ),
    gates( Qarts, Neen, Nxt, Tates ).

gates_cont( [], _Seen, [] ) :- !.
gates_cont( Parts, Seen, Gates ) :-
    gates( Parts, Seen, [], Gates ).

gates_acc( [], P, Acc, Nxt, Gates, Tates ) :-
    !,
    Gates = [P|Tates],
    Nxt = Acc.
gates_acc( ToAdd, _P, Acc, Nxt, Gates, Gates ) :-
    append( ToAdd, Acc, Nxt ).

gate_expand( [], _P, Tail, Tail ).
gate_expand( [By|T], P, Expanded, Tail ) :-
    gate_expand_pair( By, P, Expanded, ExpandedBy ),
    gate_expand( T, P, ExpandedBy, Tail ).

gate_expand_pair( Set1-Expr1, Set2-Expr2, Expanded, ExpandedBy ) :-
    % fixme: ord_disjoint/3
    ord_disjoint( Set1, Set2 ),
    !,
    ord_union( Set1, Set2, Set3 ),
    gate_expressions_norm( Expr1, Expr2, a, And ),
    gate_expressions_norm( Expr1, Expr2, o, Or  ),
    Expanded = [ Set3 - And,  % a(Expr2,Expr1),
                 Set3 - Or    % o(Expr2,Expr1)
                 | ExpandedBy ].
gate_expand_pair( _, _, Expanded, Expanded ).

gate_expressions_norm( Expr1, Expr2, Func, Expr ) :-
    gate_expression_type_args( Expr1, Func, Args1 ),
    gate_expression_type_args( Expr2, Func, Args2 ),
    append( Args1, Args2, Args ),
    sort( Args, Ord ),
    gene_expression_on_args( Ord, Func, Expr ).

gene_expression_on_args( [A,B|C], Func, Expr ) :-
    !,
    gene_expression_on_args( [B|C], Func, Right ),
    Expr =.. [Func,A,Right].
gene_expression_on_args( [A], _Func, A ).

gate_expression_type_args( Expr, Func, Args ) :-
    ( Expr =.. [Func|Args] ->
        true
        ;
        Args = [Expr]
    ).
*/

/** gate_eval( +Gn, +Gargs, +VarVecPrs, +ChildVec, -Complx, -Vector, -SubMts ).

Evaluate a boolean gate, with name Gn and args Gargs.<br>
The gate involves variables with data in VarVecPrs.<br>
Complx is the number of gates involved in Expr and the result is returned in Vector.<br>
SubMts is Pv:Ov:SubFrom;... atom giving the pairwise metrics (gate_vector_pairwise_metrics/5)<br>
for all subformulae in Gargs. ChildVec, is the target list of values.

Vects used to be a pair list of Variable-DataList pairs.

Was: gate_eval( +Expr, +VarVecPrs, -Complx, -Vec ).

*/
% gate_eval( _Set-Expr, VVs, Cx, Vec ) :- !,
    % gate_eval( Expr, VVs, Cx, Vec ).
gate_eval( Atm, [], VVs, _Against, 0, Vec, '' ) :-
    !,
    memberchk( Atm-Vec, VVs ).
gate_eval( a, Args, VVs, Against, Cx, Vec, Mts ) :- !,
    % Args = A,B
    gate_eval_args( Args, VVs, Against, Cxs, Vecs, ArgsMts ),
    % gate_eval( A, VVs, CxA, VecA ),
    % gate_eval( B, VVs, CxB, VecB ),
    % Cx is CxA + CxB + 0.99,
    % Cx is CxA + CxB + 1, % fixme: should this be length(Gars) - 1 ?
    sum_list( [1|Cxs], Cx ),
    % gate_eval_and( VecA, VecB, Vec ).
    gate_eval_vecs( Vecs, and, Vec ),
    gate_vector_pairwise_metrics( Vec, Against, Pval, Odds, Dir ),
    Gate =.. [a|Args],
    term_to_atom( Gate-m(Pval,Odds,Dir), Mt ),
    at_con( [Mt|ArgsMts], ';', Mts ).
    % gate_eval_and( Vecs, Vec ).
gate_eval( o, Args, VVs, Against, Cx, Vec, Mts ) :- !,
    gate_eval_args( Args, VVs, Against, Cxs, Vecs, ArgsMts ),
    sum_list( [1|Cxs], Cx ),
    gate_eval_vecs( Vecs, or, Vec ),
    gate_vector_pairwise_metrics( Vec, Against, Pval, Odds, Dir ),
    Gate =.. [o|Args],
    term_to_atom( Gate-m(Pval,Odds,Dir), Mt ),
    at_con( [Mt|ArgsMts], ';', Mts ).
    % gate_eval_or( Vecs, Vec ).
gate_eval( x, Args, VVs, Against, Cx, Vec, Mts ) :- !,
    gate_eval_args( Args, VVs, Against, Cxs, Vecs, ArgsMts ),
    sum_list( [1|Cxs], Cx ),
    gate_eval_vecs( Vecs, xor, Vec ),
    gate_vector_pairwise_metrics( Vec, Against, Pval, Odds, Dir ),
    Gate =.. [x|Args],
    term_to_atom( Gate-m(Pval,Odds,Dir), Mt ),
    at_con( [Mt|ArgsMts], ';', Mts ).
    % gate_eval_xor( Vecs, Vec ).
gate_eval( n, [A], VVs, Against, Cx, Vec, Mts ) :- !,
    A =.. [An|AArgs],
    gate_eval( An, AArgs, VVs, Against, CxA, VecA, ArgMt ),
    Cx is CxA + 1,
    gate_eval_not( VecA, Vec ),
    gate_vector_pairwise_metrics( Vec, Against, Pval, Odds, Dir ),
    Gate =.. [n,A],
    term_to_atom( Gate-m(Pval,Odds,Dir), Mt ),
    at_con( [Mt,ArgMt], ';', Mts ).

gate_eval_vecs( Vecs, Gate, [V|Vs] ) :-
    skim( Vecs, Vals, Rest ),
    !,
    gate_eval_values_gate( Gate, Vals, V ),
    gate_eval_vecs( Rest, Gate, Vs ).
gate_eval_vecs( _Vecs, _Gate, [] ).

gate_eval_values_gate( and, Vals, V ) :-
    gates_eval_values_and_gate( Vals, V ).
gate_eval_values_gate( or, Vals, V ) :-
    gates_eval_values_or_gate( Vals, V ).
gate_eval_values_gate( xor, Vals, V ) :-
    gates_eval_values_xor_gate( Vals, V ).

gates_eval_values_and_gate( [], 1 ).
gates_eval_values_and_gate( [0|_], 0 ) :-
    !.
gates_eval_values_and_gate( [1|T], V ) :-
    gates_eval_values_and_gate( T, V ).

gates_eval_values_or_gate( [], 0 ).
gates_eval_values_or_gate( [1|_], 1 ) :-
    !.
gates_eval_values_or_gate( [0|T], V ) :-
    gates_eval_values_or_gate( T, V ).

gates_eval_values_xor_gate( [], 0 ).
gates_eval_values_xor_gate( [1|T], Sustained ) :-
    !,
    gates_eval_values_or_gate( T, Ored ),
    gates_eval_value_not( Ored, Sustained ).
gates_eval_values_xor_gate( [0|T], V ) :-
    gates_eval_values_and_gate( T, V ).

% only used by xor...
gates_eval_value_not( 0, 1 ).
gates_eval_value_not( 1, 0 ).

/*
gate_eval( o(A,B), VVs, Cx, Vec ) :- !,
    gate_eval( A, VVs, CxA, VecA ),
    gate_eval( B, VVs, CxB, VecB ),
    Cx is CxA + CxB + 1,% fixme: should this be length(Gars) - 1 ?
    gate_eval_or( VecA, VecB, Vec ).
% experimental: 18.02.13
gate_eval( x(A,B), VVs, Cx, Vec ) :- !,
    gate_eval( A, VVs, CxA, VecA ),
    gate_eval( B, VVs, CxB, VecB ),
    Cx is CxA + CxB + 0.99,% fixme: should this be length(Gars) - 1 ?
    gate_eval_xor( VecA, VecB, Vec ).
gate_eval( n(A), VVs, Cx, Vec ) :- !,
    gate_eval( A, VVs, CxA, VecA ),
    Cx is CxA + 1,
    gate_eval_not( VecA, Vec ).
gate_eval( Atm, VVs, 0, Vec ) :-
    memberchk( Atm-Vec, VVs ).
    */

gate_eval_args( [], _VVs, _Against, [], [], [] ).
gate_eval_args( [A|As], VVs, Against, [Cx|Cxs], [V|Vs], [Mt|Mts] ) :-
    A =.. [An|AArgs],
    gate_eval( An, AArgs, VVs, Against, Cx, V, Mt ),
    gate_eval_args( As, VVs, Against, Cxs, Vs, Mts ).

/*
gate_eval_and( [], [], [] ).
gate_eval_and( [H1|T1], [H2|T2], [H|T] ) :-
    gate_eval_value_and( H1, H2, H ),
    gate_eval_and( T1, T2, T ).

gate_eval_or( [], [], [] ).
gate_eval_or( [H1|T1], [H2|T2], [H|T] ) :-
    gate_eval_value_or( H1, H2, H ),
    gate_eval_or( T1, T2, T ).

% experimental: 18.02.13
gate_eval_xor( [], [], [] ).
gate_eval_xor( [H1|T1], [H2|T2], [H|T] ) :-
    gate_eval_value_xor( H1, H2, H ),
    gate_eval_xor( T1, T2, T ).
    */

gate_eval_not( [], [] ).
gate_eval_not( [H1|T1], [H|T] ) :-
    ( H1 =:= 1 -> H = 0; H = 1 ),
    gate_eval_not( T1, T ).

gate_eval_value_and( 1, 1, 1 ) :- !.
gate_eval_value_and( _, _, 0 ).

gate_eval_value_or( 1, _, 1 ) :- !.
gate_eval_value_or( _, 1, 1 ) :- !.
gate_eval_value_or( _, _, 0 ).

gate_eval_value_xor( 1, 1, 0 ) :- !.
gate_eval_value_xor( 0, 0, 0 ) :- !.
gate_eval_value_xor( _, _, 1 ).

/* 18.02.13: i think the following is incorrect:
gates( Parts, Gates ) :-
    findall( [P]-P, member(P,Parts), PairParts ),
    gates( PairParts, [], [], Gates ).

gates( [], Seen, Acc, Gates ) :-
    gates_cont( Acc, Seen, Gates ).
gates( [P|Parts], Seen, Acc, Gates ) :-
    ( ord_memberchk(P,Seen) -> 
        Qarts = Parts,
        Neen = Seen,
        Nxt = Acc,
        Gates = Tates
        ;
        gate_expand( Parts, P, Parted, Tarted ),
        gate_expand( Seen, P, Tarted, [] ),
        gates_acc( Parted, P, Acc, Nxt, Gates, Tates ),
        P = Set-Expr,
        ( Expr = n(_) -> Parts = Qarts; Qarts = [Set-n(Expr)|Parts] ),
        ord_add_element( Seen, P, Neen )
    ),
    gates( Qarts, Neen, Nxt, Tates ).

gates_cont( [], _Seen, [] ) :- !.
gates_cont( Parts, Seen, Gates ) :-
    gates( Parts, Seen, [], Gates ).

gates_acc( [], P, Acc, Nxt, Gates, Tates ) :-
    !,
    Gates = [P|Tates],
    Nxt = Acc.
gates_acc( ToAdd, _P, Acc, Nxt, Gates, Gates ) :-
    append( ToAdd, Acc, Nxt ).

gate_expand( [], _P, Tail, Tail ).
gate_expand( [By|T], P, Expanded, Tail ) :-
    gate_expand_pair( By, P, Expanded, ExpandedBy ),
    gate_expand( T, P, ExpandedBy, Tail ).

gate_expand_pair( Set1-Expr1, Set2-Expr2, Expanded, ExpandedBy ) :-
    % fixme: ord_disjoint/3
    ord_disjoint( Set1, Set2 ),
    !,
    ord_union( Set1, Set2, Set3 ),
    gate_expressions_norm( Expr1, Expr2, a, And ),
    gate_expressions_norm( Expr1, Expr2, o, Or  ),
    Expanded = [ Set3 - And,  % a(Expr2,Expr1),
                 Set3 - Or    % o(Expr2,Expr1)
                 | ExpandedBy ].
gate_expand_pair( _, _, Expanded, Expanded ).

gate_expressions_norm( Expr1, Expr2, Func, Expr ) :-
    gate_expression_type_args( Expr1, Func, Args1 ),
    gate_expression_type_args( Expr2, Func, Args2 ),
    append( Args1, Args2, Args ),
    sort( Args, Ord ),
    gene_expression_on_args( Ord, Func, Expr ).

gene_expression_on_args( [A,B|C], Func, Expr ) :-
    !,
    gene_expression_on_args( [B|C], Func, Right ),
    Expr =.. [Func,A,Right].
gene_expression_on_args( [A], _Func, A ).

gate_expression_type_args( Expr, Func, Args ) :-
    ( Expr =.. [Func|Args] ->
        true
        ;
        Args = [Expr]
    ).
*/
