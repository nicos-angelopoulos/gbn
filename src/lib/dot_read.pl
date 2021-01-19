
% :- use_module( library(pio), [pharse_from_file/3] ). % phrase_from_file/3.
% baNaNas

%% dot_read( +File, -Graph ).
%
% An attempt to have a formal read-in grammar for graphiviz's dot language.
%
%
% @author nicos angelopoulos
% @version  0.1 2014/5/29
% @see http://www.graphviz.org/doc/info/lang.html
%
dot_read( File, Graph ) :-
	debug( dot(file), 'Reading from file: ~w', File ),
	phrase_from_file( graph(Graph), File, [] ),
    !.

graph( dot(Strict,GType,Id,Stmts) ) --> 
	dot_strict(Strict), ws, % optional
	dot_graph_type( GType ), ws,
	dot_graph_id( Id ), ws, % optional
	"{", ws, dot_stmt_list( Stmts ), ws, "}", ws.

dot_stmt_list( [Stmt|Stmts] ) --> 
	dot_stmt( Stmt ),
	{ debug(dot(statm), 'Read in statement: ~w', [Stmt] ) },
	% { (Stmt = edge([[101,120,101,99,117,116,101]|_],[]) -> trace; true) },
	dot_stmt_list_separated( Stmts ).

dot_stmt_list_optional( Stmts ) --> 
	dot_stmt_list( Stmts ).
dot_stmt_list_optional( [] ) -->  "".

dot_stmt_list_separated( Stmts ) --> 
	non_nl_ws, "\n", ws, 
	!,
	dot_stmt_list_optional( Stmts ).
dot_stmt_list_separated( Stmts ) --> 
	ws, ";", ws,
	!,
	dot_stmt_list_optional( Stmts ).
dot_stmt_list_separated( [] ) --> "".

dot_stmt( Stmt ) --> 
	dot_attr_stmt( Stmt ),
	!.
dot_stmt( Stmt ) --> 
	dot_node_stmt( Stmt ).
dot_stmt( Stmt ) --> 
	dot_edge_stmt( Stmt ).
dot_stmt( IDvar=IDval ) --> 
	dot_id( IDvar ),
	ws, "=", ws,
	dot_id( IDval ).
dot_stmt( SubG ) --> 
	dot_subgraph( SubG ).
	
dot_subgraph( subgraph(Id,Stmts) ) -->
	"subgraph", ws,
	dot_id_optional( Id ),
	"{", ws, dot_stmt_list( Stmts ), ws, "}".

dot_node_stmt( node(Id,Attrs) ) -->
	ws, dot_node_id( Id ), ws,
	dot_attr_list_optional( Attrs ),
    {debug( dot(node), 'Node: ~w', Id )}.

dot_edge_stmt( edge([Id1|OpIds],Attrs) ) -->
	ws, dot_edge_id( Id1 ), ws,
	dot_edge_rhs( OpIds ), ws,
	dot_attr_list_optional( Attrs ).

dot_attr_stmt( attr(Atype,Attrs) ) -->
	dot_attr_type( Atype ),
	% { (Atype == edge ->  trace; true ) },
	ws,
	dot_attr_list( Attrs ).

dot_attr_list_optional( Attrs ) -->
	dot_attr_list( Attrs ),
	!.
dot_attr_list_optional( [] ) --> "".

dot_attr_list( [Attrs|MAttrs] ) --> % i think the syntax allows for [xyz=1][abc=3]
	"[", ws, dot_attr_a_list_optional( Attrs ), ws, "]",
	dot_attr_list_optional( MAttrs ).

dot_attr_a_list_continuation( [Sep|TAttrs] ) --> 
	dot_attr_seperator( Sep ),
	!,
	ws,
	dot_attr_a_list( TAttrs ).
dot_attr_a_list_continuation( [] ) -->  "".

dot_attr_a_list_optional( Attrs ) --> 
	dot_attr_a_list( Attrs ),
	!.
dot_attr_a_list_optional( [] ) --> "".

dot_attr_a_list( [Idvar=Idval|TAttrs] ) -->
	dot_id( Idvar  ), ws, "=", ws, dot_id( Idval ),
	dot_attr_a_list_continuation( TAttrs ).

dot_id_optional( Id ) -->
	dot_id( Id ),
	!.
dot_id_optional( null ) --> "".

dot_id( Id ) --> 
	dot_id_string( IdCs ),
	!,
	{atom_codes( Id, IdCs )}.
dot_id( Id ) --> 
	dot_id_numeric( IdCs ),
	!,
	{ number_codes( Id, IdCs ) }.
dot_id( Id ) --> 
	dot_id_quoted( Id ),
	!.
dot_id( Id ) --> 
	dot_id_html( Id ).

dot_id_html( html(Html) ) -->
	"<",
	dot_id_html_continue( Cs ),
	{atom_codes( Html, Cs )}.

dot_id_html_continue( [] ) -->
	">",
	!.
dot_id_html_continue( [C|Cs] ) -->
	[C],
	dot_id_html_continue( Cs ).

dot_id_quoted( quote(Atom) ) -->
	"\"",
	dot_id_quoted_continue( Cs ),
	{atom_codes( Atom, Cs )}.

dot_id_quoted_continue( [] ) -->  % fixme: does not understand \"
	"\"",
	!.
dot_id_quoted_continue( [C|Cs] ) -->
	[C],
	dot_id_quoted_continue( Cs ).

dot_id_string( [C|Cs] ) --> 
	dot_id_string_starts( C ),
	dot_id_string_continues( Cs ).

dot_id_string_starts( C ) -->
	dot_id_string_character( C ),
	!.
dot_id_string_starts( C ) -->
	dot_id_underscore( C ).

dot_id_string_continues( [C|Cs] ) --> 
	dot_id_string_starts( C ),
	!,
	dot_id_string_continues( Cs ).
dot_id_string_continues( [C|Cs] ) --> 
	dot_id_digit( C ),
	!,
	dot_id_string_continues( Cs ).
dot_id_string_continues( [] ) --> "".

dot_id_string_character( C ) -->
	[C],
	{ (0'a =< C, C =< 0'z)
	  ; 
	  (0'A =< C, C =< 0'Z)
	  ; 
	  (200 =< C, C =< 377) % double check:
	  % \200-\377)
	},
	!.
dot_id_digit( C ) --> 
     [C],
	{0'0 =< C, C =< 0'9}.
% fixme: allow some options that can read HLA-A
dot_id_digit( C ) --> 
     [0'-], 
	!,
	{ C=0'- }.

dot_id_numeric( [0'-|NumCs] ) -->
	"-", 
	dot_id_number( NumCs ).
dot_id_numeric( [0'0,0'.|NumCs] ) -->
	".",
	dot_id_decimal( NumCs ).
dot_id_numeric( [C|NumCs] ) -->
	dot_id_digit( C ),
	dot_id_number( NumCs ).

dot_id_number( [0'.|Cs] ) -->
	".",
	!,
	dot_id_decimal( Cs ).
dot_id_number( [C|Cs] ) -->
	dot_id_digit( C ),
	!,
	dot_id_number( Cs ).
dot_id_number( [] ) --> "".

dot_id_decimal( [C|Cs] ) -->
	dot_id_digit( C ),
	!,
	dot_id_decimal( Cs ).
dot_id_decimal( [] ) --> "".

dot_id_underscore( 0'_ ) --> "_".

dot_edge_id( Id ) --> 
	dot_node_id( Id ).  % in the original ``subgraph'' is included

dot_edge_rhs( [Op,Id|Topid] ) -->
	dot_edge_op( Op ),
	ws, 
	dot_edge_id( Id ),
	ws,
	dot_edge_rhs_optional( Topid ).

dot_node_id( IdStr ) -->
	dot_id( Id ), 
	dot_port_optional( Id, IdStr ).

dot_port_optional( Id, Id:Port ) --> 
	dot_port( Port ),
	!.
dot_port_optional( Id, Id ) -->  "".

dot_port( Port:CompPt ) --> 
	":", dot_id( Port ),
	":", dot_compass_pt( CompPt ),
	!.
dot_port( Port ) --> 
	":", dot_id( Port ).
	
dot_edge_rhs_optional( Topid ) -->
	dot_edge_rhs( Topid ),
	!.
dot_edge_rhs_optional( [] ) --> "".

dot_strict( true ) --> "Strict".
dot_strict( false ) --> "".

dot_graph_id( ID ) -->
	dot_id_optional( ID ).

dot_edge_op( -> ) --> "->", !.
dot_edge_op( -- ) --> "--".  % checkme

dot_graph_type( graph ) --> "graph".
dot_graph_type( digraph ) --> "digraph".

dot_attr_type( graph ) --> "graph".
dot_attr_type( node ) --> "node".
dot_attr_type( edge ) --> "edge".

dot_attr_seperator( ';' ) --> ";".
dot_attr_seperator( ',' ) --> ",".

dot_compass_pt( n ) --> "n".
dot_compass_pt( ne ) --> "ne".
dot_compass_pt( e ) --> "e".
dot_compass_pt( se ) --> "se".
dot_compass_pt( s ) --> "s".
dot_compass_pt( sw ) --> "sw".
dot_compass_pt( w ) --> "w".
dot_compass_pt( nw ) --> "nw".
dot_compass_pt( c ) --> "c".
dot_compass_pt( '_' ) --> "_". % ?

non_nl_ws --> " ", !, ws.
non_nl_ws --> "\t", !, ws.
non_nl_ws --> "".

ws --> " ", !, ws.
ws --> "\t", !, ws.
ws --> "\n", !, ws.
ws --> "".
