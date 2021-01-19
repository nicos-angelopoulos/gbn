
% :- expand_file_name('$HOME/pl/lib/src',[Clean]), lib(Clean).
% :- lib( dot_read/2 ).
% :- lib( dot_read/2 ).

:- lib(os).
:- expand_file_name('$HOME/pl/lib/src',[Clean]),os_path(Clean,'graphviz/dot_read.pl',DotRead),ensure_loaded(DotRead).
:- expand_file_name('$HOME/pl/lib/src',[Clean]),os_path(Clean,'graphviz/dot_write.pl',DotRead),ensure_loaded(DotRead).

:- debug(gbn(dot_gates)).

gbn_dot_gates_defaults( [format(svg),gates(_GatesF),out(_Out),directed(false),edge_theme(jyoti)] ).
% gbn_dot_gates_defaults( [format(svg),gates(_GatesF),out(_Out),directed(false),edge_theme(four_clrs)] ).

/** gbn_dot_gates( DotFs ).
    gbn_dot_gates( +DotF, +Opts ).

    gbn_dot_gates( Opts ).

   From a BN DotF(ile) create a logic gates BN dot file, that inscribes the logic gates of the BN in DotF
   with those found in the associated GatesF (see options below).

   The output gated file is constructed from DotF by postfixing (os_postfix/3) gates to the filename.

Opts 
  * directed(Dir=false)
     whether to include edge directions

  * edge_theme(Etheme=jyoti)
     edge theme (the original can be found at four_clrs)

  * format(Fmt=svg)
     any value other than false, attempts to run upsh dot format=Fmt <lgBN.dot> on the 

  * gates(GatesF)
     when GatesF is missing, it is guessed from DotF. 
     First, a missing GatesF <DotStem>_prns_gates_best.csv and if that does not exist
     as <DotStemMinus_fclr>_prns_gates_best.csv.
     If a variable is used as input, the used file is returned.
    
  * out(Out)
     can be used to either set the output, logic-gated, dot file, or to return the created filename (var(Dot))

@author nicos angelopoulos
@version  0.1 2017/08/16
@version  0.2 2017/10/02  lib(options)

*/
gbn_dot_gates( DotF ) :-
    gbn_dot_gates( DotF, [] ).

gbn_dot_gates( DotF, Args ) :-
    options_append( gbn_dot_gates, Args, Opts ),
    options( gates(GatesF), Opts ),
    % gbn_dot_gates_csv_file( DotF, GatesF ),
    % gbn_dot_postfix_ext_file( DotF, prns_gates_best, csv, true, _, GatesF ),
    os_base_postfix_ext_old( DotF, prns_gates_best, csv, _, GatesF ),
    options( out(OutF), Opts ),
    % gbn_dot_postfix_ext_file( DotF, gated, dot, false, false, OutF ),
    os_base_postfix_ext_new( DotF, gated, dot, false, OutF ),
    debug( gbn(dot_gates), 'Picking up csv file: ~p', [GatesF] ),
    options( format(Fmt), Opts ),
    options( directed(Drc), Opts ),
    options( edge_theme(Etheme), Opts ),
    gbn_dot_gates_file( DotF, GatesF, Drc, Etheme, OutF, Fmt ).

gbn_dot_gates_file( DotF, CsvF, Drc, Etheme, OutF, Fmt ) :-
    debug_call( gbn(dot_gates), start, true ),
    dot_read( DotF, Dot ),
    Dot = dot(DotB,DotT,DotN,DotG),
    debug( gbn(dot_gates), 'Dot: ~w', Dot ),
    mtx( CsvF, Gtx, csv_read(match_arity(false)) ),
    debug_call( gbn(dot_gates), length, gates_mtx/Gtx ),
    gbn_dot_gated( Gtx, DotG, Drc, Etheme, 1, GatedG ),
    % fixme: write Gated out, and possibly display it
    Gated = dot(DotB,DotT,DotN,GatedG),
    debug( gbn(dot_gates), 'Gated: ~w', Gated ),
    % os_postfix( gated, DotF, OutF ),
    dot_write( OutF, Gated ),
    gbn_dot_disp( Fmt, OutF ),
    debug_call( gbn(dot_gates), end, true ).

gbn_dot_disp( false, _OutF ) :-
    !.
gbn_dot_disp( Fmt, OutF ) :-
    gbn_dot_disp_fmt( OutF, Fmt ).
    
gbn_dot_disp_fmt( OutF, Fmt ) :-
    absolute_file_name( path(upsh), _Path, [access(execute)] ),
    !,
    process_create( path(upsh), [dot,'format=',Fmt,OutF], [] ).
gbn_dot_disp_fmt( _OutF, _Fmt ).
    % fixme: for now just succeed, later errorr ?

os_base_postfix_ext_new( BaseF, _Psfx, _Ext, _Replace, File ) :-
    ground( File ),
    !
    os_new_from_base( File, BaseF ).
os_base_postfix_ext_new( BaseF, Psfx, Ext, Replace, File ) :-
    os_base_postfix_ext_file( BaseF, Psfx, Ext, Replace, File ),
    os_new_from_base( File, BaseF ).

os_base_postfix_ext_old( BaseF, _Psfx, _Ext, _Replace, File ) :-
    ground( File ),
    !,
    os_old_from_base( File, BaseF ).
os_base_postfix_ext_old( BaseF, Psfx, Ext, Replace, File ) :-
    os_base_postfix_ext_file( BaseF, Psfx, Ext, Replace, File ),
    exists_file( File ),
    !.
os_base_postfix_ext_old( BaseF, Psfx, Ext, Replace, File ) :-
    throw( cannot_construct_existing_from(BaseF,Psfx,Ext,Replace,File) ).

os_base_postfix_ext_file( BaseF, Psfx, Ext, Replace, File ) :-
    os_ext( _Old, Ext, BaseF, ExtF ),
    os_ext_postfix_file( Replace, ExtF, Psfx, File ).

os_ext_postfix_file( false, ExtF, Psfx, File ) :-
    os_postfix( Psfx, ExtF, File ).
os_ext_postfix_file( true, ExtF, Psfx, File ) :-
    os_postfix( Psfx, ExtF, File, replace(true) ).

os_new_from_base( File, BaseF ) :-
    exists_file( File ),
    !,
    throw( file_exists(File,from(BaseF)) ).
os_new_from_base( _File, _BaseF ).


% move somewhere in lib(options) ?
/*
gbn_dot_postfix_ext_file( _BaseF, _Psfx, _Ext, Exist, _Replace, File ) :-
    ground( File ),
    !,
    gbn_dot_postfix_ext_file_ground( Exist, error, File ).
gbn_dot_postfix_ext_file( BaseF, Psfx, Ext, Exist, Replace, File ) :-
    os_ext( _Old, Ext, BaseF, ExtF ),
    gbn_dot_postfix_file( Replace, ExtF, Psfx, File ),
    gbn_dot_postfix_ext_file_ground( Exist, fail, File ),
    !.
gbn_dot_postfix_ext_file( BaseF, _Psfx, _Ext, _Exist, _Replace, _File ) :-
    throw( cannot_locate_gates_file_for(BaseF) ).

gbn_dot_postfix_file( false, ExtF, Psfx, File ) :-
    os_postfix( Psfx, ExtF, File ).
gbn_dot_postfix_file( true, ExtF, Psfx, File ) :-
    os_postfix( Psfx, ExtF, File, replace(true) ).
*/

/*
gbn_dot_gates_csv_file( DotF, CsvF ) :-
    ground( CsvF ),
    exists_file( CsvF ).
gbn_dot_gates_csv_file( DotF, CsvF ) :-
    os_ext( dot, csv, DotF, CsvShF ),
    os_postfix( prns_gates_best, CsvShF, CsvF ),
    exists_file( CsvF ),
    !.
gbn_dot_gates_csv_file( DotF, CsvF ) :-
    os_ext( dot, csv, DotF, CsvShF ),
    os_postfix( prns_gates_best, CsvShF, CsvF, replace(true) ),
    exists_file( CsvF ),
    !.
gbn_dot_gates_csv_file( DotF, _CsvF ) :-
    throw( cannot_constuct_existing_csv_file_name_from_dot_file(DotF) ).
    */

/*
gbn_dot_postfix_ext_file_ground( true, AtFail, ExistsF ) :- 
    gbn_dot_postfix_ext_file_ground_exists( ExistsF, AtFail ).
gbn_dot_postfix_ext_file_ground( false, AtFail, ExistsF ) :- 
    gbn_dot_postfix_ext_file_ground_does_not_exist( ExistsF, AtFail ).

gbn_dot_postfix_ext_file_ground_does_not_exist( ExistsF, AtFail ) :- 
    exists_file( ExistsF ),
    !,
    gbn_dot_postfix_ext_file_ground_does_not_exist_fail( AtFail, ExistsF ).
gbn_dot_postfix_ext_file_ground_does_not_exist( _ExistsF, _AtFail ).

gbn_dot_postfix_ext_file_ground_does_not_exist_fail( error, ExistsF ) :-
    throw( ground_output_file_does_not_exist(ExistsF) ).
% gbn_dot_postfix_ext_file_ground_does_not_exist_fail( fail, ExistsF ) :- fail.

gbn_dot_postfix_ext_file_ground_exists( ExistsF, _AtFail ) :- 
    exists_file( ExistsF ),
    !.
gbn_dot_postfix_ext_file_ground_exists( ExistsF, AtFail ) :- 
    gbn_dot_postfix_ext_file_ground_exists_fail( AtFail, ExistsF ).

gbn_dot_postfix_ext_file_ground_exists_fail( error, ExistsF ) :-
    throw( output_file_already_exists(ExistsF) ).
% gbn_dot_postfix_ext_file_ground_exists_fail( fail, _ExistsF ) :- fail. % PS !
*/

/** gbn_dot_gated( +Gtx, +Dot, +Drc, +NxtLG, -Gated ).

Gated 

*/
gbn_dot_gated( [], Gated, _Drc, _Etheme, _, Gated ).
gbn_dot_gated( [VarsRow,GateRow|T], Dot, Drc, Etheme, CntLG, Gated ) :-
    VarsRow =.. [_Rnm|Vars],
    % GateRow = row(_,_,_Dir,_Pv,_LgR,GateAtm),
    arg( 6, GateRow, GateAtm ),
    term_to_atom( Gate, GateAtm ),
    debug( gbn(dot_gates), 'Vars: ~w. Gate: ~w', [Vars,Gate] ),
    gbn_dot_gate( Vars, Gate, GateRow, Dot, Drc, Etheme, CntLG, MidDot, NxtLG ),
    gbn_dot_gated( T, MidDot, Drc, Etheme, NxtLG, Gated ).

gbn_dot_gate( [Child|Pas], Gate, GateRow, Graph, Drc, Etheme, CntLG, NxtGraph, NxtLG ) :-
    gbn_dot_edges_sel( Pas, Child, Graph, Drc, RedGraph, SelEs ),
    debug( gbn(dot_gates), 'Child: ~w. Parents: ~w. Sel.Edges: ~w', [Child,Pas,SelEs] ),
    gbn_dot_gate_add( Gate, GateRow, Pas, Child, RedGraph, Drc, Etheme, SelEs, CntLG, NxtGraph, Rangles, NxtLG ),
    ( Rangles == [] -> true; throw( remaing_dangles(Rangles) ) ).

gbn_dot_gate_add( n(OpN), GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Langles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gate_add( OpN, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, NxtLG ),
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(circle), % fixme: circle
                               fillcolor=quote('#87A96B'), fixedsize=true, width=0.2, height=0.2
                              ]] ),
    % fixedsize=true, width=0.3, height=0.3]
    gbn_dot_gate_add_node( LGraph, Node, NxtGraph ).
gbn_dot_gate_add( o(OpL,OpR), GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gate_add( OpL, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, LLG ),
    gbn_dot_gate_add( OpR, row(0,0,none,0.9,1), Pas, Gid, LGraph, Drc, Etheme, Langles, LLG, RGraph, Rangles, NxtLG ),
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(invhouse),
                               fillcolor=quote('#87A96B'), fixedzise=true, width=0.4, height=0.3
                              ]] ),
    % fixedsize=true, width=0.3, height=0.3]
    gbn_dot_gate_add_node( RGraph, Node, NxtGraph ).
gbn_dot_gate_add( a(OpL,OpR), GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gate_add( OpL, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, LLG ),
    gbn_dot_gate_add( OpR, row(0,0,none,0.9,1), Pas, Gid, LGraph, Drc, Etheme, Langles, LLG, RGraph, Rangles, NxtLG ),
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(invtriangle),
                               fillcolor=quote('#87A96B'), fixedsize=true, width=0.4, height=0.4
                              ]] ),
    % fixedsize=true, width=0.3, height=0.3]
    gbn_dot_gate_add_node( RGraph, Node, NxtGraph ).
% gbn_dot_gate_add( a(OpL,OpR), GateRow, Pas, Child, Graph, Dangles, CntLG, NxtGraph, NxtLG ) :-
gbn_dot_gate_add( x(OpL,OpR), GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gate_add( OpL, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, LLG ),
    gbn_dot_gate_add( OpR, row(0,0,none,0.9,1), Pas, Gid, LGraph, Drc, Etheme, Langles, LLG, RGraph, Rangles, NxtLG ),
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(invtrapezium),
                               fillcolor=quote('#87A96B'), fixedzise=true, width=0.4, height=0.3
                              ]] ),
    % fixedsize=true, width=0.3, height=0.3]
    gbn_dot_gate_add_node( RGraph, Node, NxtGraph ).

% gbn_dot_gate_add( n(Op), GateRow, Pas, Child, Graph, Dangles, CntLG, NxtGraph, NxtLG ) :-

gbn_dot_gate_add( Pa, _GateRow, Pas, Child, Graph, Drc, _Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-
    % by not using Etheme, we let the graph inherit the properties of its ancestor
    memberchk( Pa, Pas ),
    !,
    ( ( select( edge([quote(Pa),->,quote(Child)],Atts),Dangles,Rangles)
        ;
        select( edge([quote(Pa),--,quote(Child)],Atts),Dangles,Rangles)
        )  ->
        true
        ;
        throw( not_in_dangles(Pa,Dangles) )
    ),
    NxtLG = CntLG,
    gbn_direction_to_edge_atom( Drc, Arrow ),
    gbn_dot_gate_add_edge( Graph, edge([quote(Pa),Arrow,quote(Child)],Atts), NxtGraph ).
gbn_dot_gate_add( Other, GateRow, Pas, Child, _Graph, _Drc, _Etheme, _Dangles, _CntLG, _NxtGraph, _NxtLG ) :-
    throw( cannot_add_gate(Other,Pas,Child,GateRow) ).

gbn_dot_gate_add_from_gate_edge( From, To, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ) :-
    atom_concat( gbn_dot_gate_dir_pvalue_edge_atts_, Etheme, Pname ),
    Goal =.. [Pname,GDir,GPvl,Ovl,EdgeAtts],
    call( Goal ),
    gbn_direction_to_edge_atom( Drc, Arrow ),
    Edge = edge([quote(From),Arrow,quote(To)],EdgeAtts),
    gbn_dot_gate_add_edge( Graph, Edge, PreGraph ).

gbn_direction_to_edge_atom( false, '--' ).
gbn_direction_to_edge_atom( true, '->' ).

gbn_dot_gate_dir_pvalue_edge_atts_jyoti( GDir, Pvl, Ovl, Atts ) :-
    % gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( GDir, Pvl, EdgeClr ),
    ( (Pvl < 0.05;GDir==none) -> PenStyle=solid; PenStyle=dashed ),
    gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( GDir, 0.0001, EdgeClr ),
    ( (GDir == none;GDir== null) ->
        Pw is 1
        ;
        ( Ovl == '1.0Inf' -> 
            Pw = 4
            ;
            ( Ovl > 1 ->
                Pw is 1 + min(3,log(min(2*Ovl,10)))
                ;
                Pw is 1 + min(3,exp(0.6-Ovl))
            )
        )
    ),
    Atts = [[style=quote(PenStyle),color=quote(EdgeClr),penwidth(Pw)]].

gbn_dot_gate_dir_pvalue_edge_atts_four_clrs( GDir, Pvl, _Ovl, Atts ) :-
    gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( GDir, Pvl, EdgeClr ),
    Atts = [[color=quote(EdgeClr)]].

gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( null, _Pvl, black ).
gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( none, _Pvl, black ).
gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( mut_excl, Pvl, Clr ) :-
    Signif = '#BF812D',
    NonSignif = '#CC0000',
    gbn_dot_gate_pvalue_edge_colour( Pvl, Signif, NonSignif, Clr ).
gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( co_occur, Pvl, Clr ) :-
    Signif = '#35978F',
    NonSignif = '#4B5320',
    gbn_dot_gate_pvalue_edge_colour( Pvl, Signif, NonSignif, Clr ).

gbn_dot_gate_pvalue_edge_colour( Pvl, Signif, _NonSignif, Clr ) :-
    Pvl < 0.05,
    !,
    Signif = Clr.
gbn_dot_gate_pvalue_edge_colour( _Pvl, _Signif, NonSignif, Clr ) :-
    NonSignif = Clr.

/** gbn_dot_gate_add_edge( +Graph, +Edge, -EGraph ).

    Add Edge to Graph to produce EGraph. The edge is inserted just before 
    the first existing edge/2 term in Graph (list) or as the last item in Graph.
    
*/
gbn_dot_gate_add_edge( [], Edge, [Edge] ).
gbn_dot_gate_add_edge( [edge(Link,Atts)|T], Edge, [Edge,edge(Link,Atts)|T] ) :-
    !.
gbn_dot_gate_add_edge( [H|T], Edge, [H|M] ) :-
    gbn_dot_gate_add_edge( T, Edge, M ).

gbn_dot_gate_add_node( [], Node, [Node] ).
gbn_dot_gate_add_node( [node(A,B)|T], Node, [Node,node(A,B)|T] ) :-
    !.
gbn_dot_gate_add_node( [H|T], Node, [H|R] ) :-
    gbn_dot_gate_add_node( T, Node, R ).

/** gbn_dot_gate_count( +CntLG, -Gid, -NxtLG ).

    Construct the next valid gate id, and increase the gate counter.

*/
gbn_dot_gate_count( CntLG, Gid, NxtLG ) :-
    NxtLG is CntLG + 1,
    atomic_list_concat( [gate,CntLG], '_', Gid ).

/** gbn_dot_edges_sel( +Pas, +Child, +Graph, +Drc, -GraphNE, -SelEs ).


Select each directed edge in Graph from parents in Pas list to node Child
producing GraphNE, the reduced graph and the selected edges with Child replaced by an
unbound variable.

*/
gbn_dot_edges_sel( [], _Child, Graph, _Drc, Graph, [] ).
gbn_dot_edges_sel( [Pa|Pas], Child, Graph, Drc, NxtGraph, [SelE|TSelEs] ) :-
    gbn_direction_to_edge_atom( Drc, Arrow ),
    gbn_dot_edge_sel( Pa, Child, Graph, Arrow, GraphNE, SelE ),
    gbn_dot_edges_sel( Pas, Child, GraphNE, Drc, NxtGraph, TSelEs ).

gbn_dot_edge_sel( Pa, Child, Graph, Arrow, GraphNE, SelE ) :-
    select( edge([Pa,->,Child],Atts), Graph, GraphNE ),
    !,
    SelE = edge([quote(Pa),Arrow,_],Atts).
gbn_dot_edge_sel( Pa, Child, Graph, Arrow, GraphNE, SelE ) :-
    select( edge([quote(Pa),->,quote(Child)],Atts), Graph, GraphNE ),
    !,
    SelE = edge([quote(Pa),Arrow,_],Atts).
gbn_dot_edge_sel( Pa, Child, Graph, Arrow, GraphNE, SelE ) :-
    select( edge([Pa,--,Child],Atts), Graph, GraphNE ),
    !,
    SelE = edge([quote(Pa),Arrow,_],Atts).
gbn_dot_edge_sel( Pa, Child, Graph, Arrow, GraphNE, SelE ) :-
    select( edge([quote(Pa),--,quote(Child)],Atts), Graph, GraphNE ),
    !,
    SelE = edge([quote(Pa),Arrow,_],Atts).
gbn_dot_edge_sel( Pa, Child, Graph, _Arrow, _GraphNE, _SelE ) :-
    throw( cannot_find_edge_fron_to_in(Pa,Child,Graph) ).
