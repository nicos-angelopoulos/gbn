
% :- expand_file_name('$HOME/pl/lib/src',[Clean]), lib(Clean).
% :- lib(dot_read/2).
% :- lib(dot_read/2).

:- lib(os_lib).
:- lib(dot/1).
:- lib(dot_read/2).
:- lib(dot_write/2).

gbn_gates_net_defaults( Defs ) :-
        Defs = [    format(svg), gates(_GatesF),
                    out(_Out), directed(false), 
                    edge_theme(jyoti) 
              ].

/** gbn_gates_net( DotF ).
    gbn_gates_net( +DotF, +Opts ).

   From a BN DotF(ile) create a logic gates BN dot file, that inscribes the logic gates of the BN in DotF
   with those found in the associated GatesF (see options below). All fitted gates are already in the per-family
   gate files, GatesF, here we simply chose the "best" gate from each family.

   Currently, this is the gate with smallest p.value. Ignoring direction seems the best default.

   If no output file is given, it is constructed from DotF by postfixing (os_postfix/3) gates_best to the filename.

Opts 
  * directed(Dir=false)
     whether to include edge directions

  * edge_theme(Etheme=jyoti)
     edge theme (the original can be found at four_clrs)

  * format(Fmt=svg)
     any value other than false, attempts to run upsh dot format=Fmt <lgBN.dot> on the generated dot file

  * gates(GatesF)
     when GatesF is missing, it is guessed from DotF. 
     First, a missing GatesF <DotStem>_gates_best.csv and if that does not exist
     as <DotStemMinus_fclr>_gates_best.csv.
     If a variable is used as input, the used file is returned.
    
  * out(Out)
     can be used to either set the output, logic-gated, dot file, or to return the created filename (var(Dot))

@author nicos angelopoulos
@version  0.1 2017/08/16
@version  0.2 2017/10/02  lib(options)

*/
gbn_gates_net( DotF ) :-
    gbn_gates_net( DotF, [] ).

gbn_gates_net( DotF, Args ) :-
    options_append( gbn_gates_net, Args, Opts ),
    options( gates(GatesF), Opts ),
    % gbn_dot_gates_csv_file( DotF, GatesF ),
    % gbn_dot_postfix_ext_file( DotF, prns_gates_best, csv, true, _, GatesF ),
    os_base_postfix_ext_old( DotF, gates_best, csv, _, GatesF ),
    options( out(OutF), Opts ),
    % gbn_dot_postfix_ext_file( DotF, gated, dot, false, false, OutF ),
    os_base_postfix_ext_new( DotF, gated, dot, false, OutF ),
    debug( gbn(gates_net), 'Picking up csv file: ~p', [GatesF] ),
    % ( GatesF == 'dutch_driver_muts_min15-18.03.01/dutch_driver_muts_min15-e1_gates_best.csv' -> spy(gbn_dot_gate_add) ; true ),
    options( format(Fmt), Opts ),
    options( directed(Drc), Opts ),
    options( edge_theme(Etheme), Opts ),
    gbn_gates_net_file( DotF, GatesF, Drc, Etheme, OutF, Fmt ).

gbn_gates_net_file( DotF, CsvF, Drc, Etheme, OutF, Fmt ) :-
    dot_read( DotF, Dot ),
    Dot = dot(DotB,DotT,DotN,DotG),
    % debug( gbn(gates_net), 'Dot: ~w', Dot ),
    % mtx( CsvF, Gtx, csv_read(match_arity(false)) ),
    mtx( CsvF, Gtx, [match(false),convert(true)] ),
    gbn_dot_gated( Gtx, DotG, Drc, Etheme, 1, GatedG ),
    % fixme: write Gated out, and possibly display it
    Gated = dot(DotB,DotT,DotN,GatedG),
    % debug( gbn(gates_net), 'Gated: ~w', Gated ),
    % os_postfix( gated, DotF, OutF ),
    dot_write( OutF, Gated ),
    gbn_dot_disp( Fmt, OutF ).

gbn_dot_disp( false, _OutF ) :-
    !.
gbn_dot_disp( Fmt, OutF ) :-
    gbn_dot_disp_fmt( OutF, Fmt ).
    
gbn_dot_disp_fmt( OutF, Fmt ) :-
    % absolute_file_name( path(upsh), _Path, [access(execute)] ),
    % !,
    % process_create( path(upsh), [dot,'format=',Fmt,OutF], [] ).
    dot_file( Fmt, OutF ),
    !.
gbn_dot_disp_fmt( _OutF, _Fmt ).
    % fixme: for now just succeed, later error ?

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


/** gbn_dot_gated( +Gtx, +Dot, +Drc, +NxtLG, -Gated ).

Gated 

*/
gbn_dot_gated( [], Gated, _Drc, _Etheme, _, Gated ).
gbn_dot_gated( [VarsRow,GateRow|T], Dot, Drc, Etheme, CntLG, Gated ) :-
    VarsRow =.. [_Rnm|Vars],
    % GateRow = row(_,_,_Dir,_Pv,_LgR,GateAtm),
    arg( 6, GateRow, GateAtm ),
    term_to_atom( Gate, GateAtm ),
    debug( gbn(gates_net), 'Vars: ~w. Gate: ~w', [Vars,Gate] ),
    debug( gbn(gates_net), 'GateAtm: ~w', [GateAtm] ),
    gbn_dot_gate( Vars, Gate, GateRow, Dot, Drc, Etheme, CntLG, MidDot, NxtLG ),
    gbn_dot_gated( T, MidDot, Drc, Etheme, NxtLG, Gated ).

gbn_dot_gate( [Child|Pas], Gate, GateRow, Graph, Drc, Etheme, CntLG, NxtGraph, NxtLG ) :-
    gbn_dot_edges_sel( Pas, Child, Graph, Drc, RedGraph, SelEs ),
    debug( gbn(gates_net), 'Child: ~w. Parents: ~w. Sel.Edges: ~w', [Child,Pas,SelEs] ),
    Gate =.. [Gn|GArgs],
    gbn_dot_gate_add( Gn, GArgs, GateRow, Pas, Child, RedGraph, Drc, Etheme, SelEs, CntLG, NxtGraph, Rangles, NxtLG ),
    ( Rangles == [] -> true; throw( remaining_dangles(Rangles) ) ).

gbn_gate_dot_node( Gn, Gid, Node ) :-
    gbn_gate_dot_atts( Gn, Atts ),
    Node = node( quote(Gid), Atts ).

gbn_gate_dot_atts(a, [[style=filled,label=quote(''),shape=quote(invtriangle),fillcolor=quote('#87A96B'),fixedsize=true,width=0.4,height=0.4]] ).
gbn_gate_dot_atts(o, [[style=filled,label=quote(''),shape=quote(invhouse),fillcolor=quote('#87A96B'),fixedzise=true,width=0.4,height=0.3]] ).
gbn_gate_dot_atts(n, [[style=filled,label=quote(''),shape=quote(circle),fillcolor=quote('#87A96B'),fixedsize=true,width=0.2,height=0.2]] ).
gbn_gate_dot_atts(x, [[style=filled,label=quote(''),shape=quote(invtrapezium),fillcolor=quote('#87A96B'),fixedzise=true,width=0.4,height=0.3]] ).

/* 18.03.07: commented out as it should be dealt with the clause below...also it never is the case that n() will be the outer layer
gbn_dot_gate_add( n, [OpN], GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Langles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    % here: deconstruct arg( 7 which is the new Gate-m(Pv,Ov,Dir)
        arg( 7, GateRow, SubGateMtsAtm ),
        % term_to_atom( SubGateMts, SubGateMtsAtm ),
        % memberchk( OpN-m(OpN) )

    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    OpN =.. [OpNNm|OpNArgs],
    gbn_gate_direc_complement( GDir, GDirCmp ),
    Mow = row(0,0,GDirCmp,GPvl,Ovl,0,SubGateMtsAtm),
    gbn_dot_gate_add( OpNNm, OpNArgs, Mow, Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, NxtLG ),
    gbn_gate_dot_node( n, Gid, Node ),
    % fixedsize=true, width=0.3, height=0.3]
    gbn_dot_gate_add_node( LGraph, Node, NxtGraph ).
    */
% gbn_dot_gate_add( o(OpL,OpR), GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-

gbn_dot_gate_add( Gn, [Opea|Opeas], GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Oangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    arg( 7, GateRow, SubGateMtsAtm ),
    % term_to_atom( SubGateMtsTerm, SubGateMtsAtm ),
    atomic_list_concat( SubGateMtsAtms, ';', SubGateMtsAtm ),
    maplist( term_to_atom, SubGateMts, SubGateMtsAtms ),

    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gates_add( [Opea|Opeas], Pas, Gid, SubGateMts, SubGateMtsAtm, PreGraph, Drc, Etheme, Dangles, PreLG, OGraph, Oangles, NxtLG ),
    %
    % gbn_dot_gate_add( OpL, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, LLG ),
    % gbn_dot_gate_add( OpR, row(0,0,none,0.9,1), Pas, Gid, LGraph, Drc, Etheme, Langles, LLG, RGraph, Rangles, NxtLG ),
    %
    gbn_gate_dot_node( Gn, Gid, Node ),
    gbn_dot_gate_add_node( OGraph, Node, NxtGraph ).

/*
gbn_dot_gate_add( a, Opeas, GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Oangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gates_add( Opeas, Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, OGraph, Oangles, NxtLG ),
    gbn_gate_dot_shape( Gn, Node ),
    gbn_dot_gate_add_node( OGraph, Node, NxtGraph ).

gbn_dot_gate_add( x, Opeas, GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Oangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gates_add( Opeas, Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, OGraph, Oangles, NxtLG ),
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(invtrapezium), fillcolor=quote('#87A96B'), fixedzise=true, width=0.4, height=0.3 ]] ),
    gbn_dot_gate_add_node( OGraph, Node, NxtGraph ).
*/

    /*
gbn_dot_gate_add( a(OpL,OpR), GateRow, Pas, Child, Graph, Drc, Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-
    !,
    gbn_dot_gate_count( CntLG, Gid, PreLG ),
    arg( 3, GateRow, GDir ),
    arg( 4, GateRow, GPvl ),
    arg( 5, GateRow, Ovl ),
    gbn_dot_gate_add_from_gate_edge( Gid, Child, GDir, GPvl, Ovl, Graph, Drc, Etheme, PreGraph ),
    gbn_dot_gate_add( OpL, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, LLG ),
    gbn_dot_gate_add( OpR, row(0,0,none,0.9,1), Pas, Gid, LGraph, Drc, Etheme, Langles, LLG, RGraph, Rangles, NxtLG ),
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(invtriangle), fillcolor=quote('#87A96B'), fixedsize=true, width=0.4, height=0.4 ]] ),
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
    Node = node( quote(Gid), [[style=filled, label=quote(''),shape=quote(invtrapezium), fillcolor=quote('#87A96B'), fixedzise=true, width=0.4, height=0.3 ]] ),
    % fixedsize=true, width=0.3, height=0.3]
    gbn_dot_gate_add_node( RGraph, Node, NxtGraph ).
    */

% gbn_dot_gate_add( n(Op), GateRow, Pas, Child, Graph, Dangles, CntLG, NxtGraph, NxtLG ) :-
gbn_dot_gate_add( Pa, [], _GateRow, Pas, Child, Graph, Drc, _Etheme, Dangles, CntLG, NxtGraph, Rangles, NxtLG ) :-
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
gbn_dot_gate_add( Other, OArgs, GateRow, Pas, Child, _Graph, _Drc, _Etheme, _Dangles, _CntLG, _NxtGraph, _NxtLG ) :-
    throw( cannot_add_gate(Other,OArgs,Pas,Child,GateRow) ).

gbn_dot_gates_add( [], _Pas, _Gid, _SubGateMts, _SubGateMtsAtm, Graph, _Drc, _Etheme, Dangles, LG, Graph, Dangles, LG ).
gbn_dot_gates_add( [Opea|Opeas], Pas, Gid, SubGateMts, SubGateMtsAtm, PreGraph, Drc, Etheme, PreDangles, PreLG, Graph, Dangles, LG ) :-
    % gbn_dot_gate_add( OpL, row(0,0,none,0.9,1), Pas, Gid, PreGraph, Drc, Etheme, Dangles, PreLG, LGraph, Langles, LLG ),
    Opea =.. [Gn|GArgs],
    ( memberchk(Opea-m(OpPv,OpOv,OpDrc),SubGateMts) ->
        GateRow = row(0,0,OpDrc,OpPv,OpOv,null_expr,SubGateMtsAtm)
        ;
        GateRow = row(0,0,none,0.9,1,null_expr,SubGateMtsAtm)
    ),
    gbn_dot_gate_add( Gn, GArgs, GateRow, Pas, Gid, PreGraph, Drc, Etheme, PreDangles, PreLG, NxtGraph, NxtDangles, NxtLG ),
    gbn_dot_gates_add( Opeas, Pas, Gid, SubGateMts, SubGateMtsAtm, NxtGraph, Drc, Etheme, NxtDangles, NxtLG, Graph, Dangles, LG ).

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
    % ( (Pvl < 0.05;GDir==none) -> PenStyle=solid; PenStyle=dashed ),
    ( Pvl < 0.05 -> PenStyle=solid; PenStyle=dashed ),
    gbn_dot_gate_dir_pvalue_edge_atts_four_clrs_clr( GDir, 0.0001, EdgeClr ),
    ( (GDir == none;GDir== null) ->
        Pw is 1
        ;
        gbn_fisher_edge_width( prop, Ovl, Pw )
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

gbn_gate_direc_complement( co_occur, mut_excl ) :- !.
gbn_gate_direc_complement( mut_excl, co_occur ) :- !.
gbn_gate_direc_complement( A, A ).
