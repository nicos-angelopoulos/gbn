%% dot_write( +File, +Graph ).
%
% Graph should be as those produced by read_dot/2.
%
%
% @author nicos angelopoulos
% @version  0.1 2014/5/29
%
dot_write( File, Graph ) :-
	phrase( dot_spit(Graph), String ),
	setup_call_cleanup(
	    open(File, write, Out ),
	    format(Out, '~s', [String]),
	    close(Out) ).

dot_spit( dot(Strict,Gtype,Id,Stmts) ) -->
	dot_spit_graph_strict( Strict ), 
	dot_spit_atom( Gtype ), " ",
	dot_spit_atom( Id ),
	" {", "\n",
	dot_spit_statements( Stmts ),
	"}".

dot_spit_statements( [] ) --> "".
dot_spit_statements( [H|T] ) -->
	"\t",
	{ debuc(dot_write, 'Statement:~w', H )},
	dot_spit_statement( H ),
	";\n",
	dot_spit_statements( T ).

dot_spit_statement( attr(Atype,Attrs) ) -->
	dot_spit_atom( Atype ), 
	" ",
	dot_spit_attrs( Attrs ).

dot_spit_statement( Var=Val ) -->
	dot_spit_id( Var ),
	"=",
	dot_spit_id( Val ).

dot_spit_statement( node(Id,Attrs) ) -->
	dot_spit_id( Id ),
	dot_spit_attrs( Attrs ).
dot_spit_statement( edge(IdOpChain,Attrs) ) -->
	dot_spit_id_op_chain( IdOpChain ),
	dot_spit_attrs( Attrs ).
dot_spit_statement( subgraph(Id,Stmts) ) -->
	"subgraph ",
	dot_spit_id( Id ),
	"{",
	dot_spit_statements( Stmts ),
	"}".
dot_spit_statement( nl ) --> "\n".

dot_spit_id_op_chain( [Id|T] ) --> 
	dot_spit_id( Id ),
	dot_spit_id_op_chain_rhs( T ).

dot_spit_id_op_chain_rhs( [] ) --> "".
dot_spit_id_op_chain_rhs( [Op,Id|T] ) --> 
	dot_spit_edge_op( Op ),
	dot_spit_id( Id ),
	dot_spit_id_op_chain_rhs( T ).

dot_spit_id( quote(Id) ) --> 
	 { !, atom_codes(Id,IdCs) },
	 "\"",IdCs, "\"".
dot_spit_id( html(Id) ) --> 
	 { !, atom_codes(Id,IdCs) },
	 "<",IdCs, ">".
dot_spit_id( Id ) --> 
	 {atom(Id)},
	 !,
	 dot_spit_id_atom( Id ).
dot_spit_id( Id ) --> 
	{number(Id),
	 number_codes(Id,IdCs) },
	IdCs.

% fixme: more specials? 
dot_spit_id_atom( Id ) -->
	{Specials = ['-'], member(Spc,Specials), sub_atom(Id,_,_,_,Spc),
	 atom_codes(Id,IdCs) },
	!,
	[0'"],
	IdCs,
	[0'"].
dot_spit_id_atom( Id ) -->
	{atom_codes(Id,IdCs) },
	IdCs.

dot_spit_attrs( InAttrs ) -->
	{to_nest_list( InAttrs, Attrs )},
	dot_spit_nest_attrs( Attrs ).

dot_spit_nest_attrs( [] ) --> "".
dot_spit_nest_attrs( [Al|Als] ) -->
	"[",
	dot_spit_attrs_list( Al ),
	"]",
	dot_spit_nest_attrs( Als ).

dot_spit_attrs_list( [] ) --> "".
dot_spit_attrs_list( [A|As] ) -->
	dot_spit_attr( A ),
	dot_spit_attrs_sep( As, Bs ),
	dot_spit_attrs_list( Bs ).

dot_spit_attr( Var=Val ) --> 
    !,
	dot_spit_id( Var, Val ).
dot_spit_attr( Comp ) -->      % be relaxed
	{ compound( Comp ),
	   ground( Comp ),
	   Comp =.. [Var,Val] },
	dot_spit_id( Var, Val ).

dot_spit_id( Var, Val ) -->
	dot_spit_id( Var ),
	"=",
	dot_spit_id( Val ).

dot_spit_attrs_sep( [Sep|Bs], Bs ) -->
	dot_spit_attrs_sep_token( Sep ),
	" ",
    !.
% dot_spit_attrs_sep( Bs, Bs ) -->    "".     % be relaxed about missing seps
dot_spit_attrs_sep( [], [] ) -->  !, "".
dot_spit_attrs_sep( Bs, Bs ) -->  ",".

dot_spit_attrs_sep_token( ',' ) --> ",".
dot_spit_attrs_sep_token( ';' ) --> ";".

dot_spit_graph_strict( false ) --> [].
dot_spit_graph_strict( true ) --> ["strict "].

dot_spit_edge_op( -> ) --> " -> ".
dot_spit_edge_op( -- ) --> " -- ".

dot_spit_atom( Gtype ) --> 
	{atom_codes(Gtype,GtypeCs)},
	GtypeCs.

to_nest_list( [H|T], Nest ) :-
	is_list(H),
	!,
	Nest = [H|T].
to_nest_list( Flat, [Flat] ).
