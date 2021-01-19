
% :- lib(real).
% :- lib(options).
% :- lib(disp_bn).
% :- lib(stoics_lib:kv_decompose/3).

:- lib(cross_column_fisher_test/6).
:- lib(gbn_mtx_df/6).
:- lib(r("RColorBrewer")).

% gbn_fisher_theme_colours( +Theme, -Cco, -Ccs, -Cex, -Ces ).
% 
% Colours for Fisher edges.
%
% Cco - Co-occurance (non significant)
%
% Ccs - Co-occurance (significant)
% 
% Cex - Mutual exclusivity (non significant)
%
% Ces - Mutual exclusivity (significant)
%
gbn_fisher_theme_colours( org, "#4B5320", "#35978F", "#CC0000", "#BF812D" ).
gbn_fisher_theme_colours( bic, "#35978F", "#35978F", "#BF812D", "#BF812D" ).

gbn_fisher_theme_dashed( org, solid ).
gbn_fisher_theme_dashed( bic, dashed ).

gbn_fisher_theme_edge_width( org, fix ).
gbn_fisher_theme_edge_width( bic, prop ).

gbn_fisher_theme_node_pen_width( org, fix ).
gbn_fisher_theme_node_pen_width( bic, rel_prop ).

% gbn_fisher_theme_background( org, '#D3D3D3' ).
gbn_fisher_theme_background( org, '#FFFFFF' ).
gbn_fisher_theme_background( bic, '#FFFFFF' ).

% fixme: add rvar ? 
gbn_fisher_net_defaults( Args, Defs ) :- 
    DefTheme = bic,
    Defs    = [ 
            directed(false),
            % bground('#D3D3D3'),
            bground(Bg),
            postfix(Psfx),
            adjust('BH'),
            clr_theme(DefTheme),
            clr_co(Cco),
            clr_co_signif(Ccs),
            clr_ex(Cex),
            clr_ex_signif(Ces),
            dashed(Dashed),
            edge_width(Ew),
            node_pen_width(Nw),
            edge_metrics(_)   % return value option
        ],
    ( memberchk(clr_theme(Theme),Args) -> true; Theme = DefTheme ),
    ( Theme == DefTheme -> Psfx = ''; Psfx = Theme ),
    atomic( Theme ),
    % fixme: use known/n predicate
    gbn_fisher_theme_dashed( Theme, Dashed ),
    gbn_fisher_theme_edge_width( Theme, Ew ),
    gbn_fisher_theme_node_pen_width( Theme, Nw ),
    gbn_fisher_theme_background( Theme, Bg ),
    gbn_fisher_theme_colours( Theme, Cco, Ccs, Cex, Ces ).

/** gbn_fisher_net( +DatF, +BnF, -FisheredF ).
    gbn_fisher_net( +DatF, +BnF, -FisheredF, +Opts ).

Colour  Gobnilp BNs with co-occurance / mutual exclusivity scores
according to Fisher test. 

Options Opts are also passed to disp_bn/2.

Opts
  * adjust(Adj='BH')
    either false, or a method recognised by R's p.adjust() function

  * bground(ClrBg="#D3D3D3")
    background colour

  * clr_co(ClrCo)
    colour for co-occurance (non significant level), default depends on Theme

  * clr_co_signif(ClrCoSgn)
    colour for co-occurance (for significant values), default depends on Theme

  * clr_ex(ClrEx)
    colour for mutual exclusivity (non significant level), default depends on Theme

  * clr_ex_signif(ClrExSgn)
    colour for mutual exclusivity (for significant values), default depends on Theme

  * clr_theme(Theme=bic)
    defines default colour values.
    =org= theme is:
     * ClrCo="#4B5320" 
        (green)
     * ClrCoSgn="#35978F" 
        (bluegreen)
     * ClrEx="#CC0000" 
        (red)
     * ClrExSgn="#BF812D" 
        (golden)
     * Dash=solid 
       and
     * Ew=fix

    =bic= theme is:
     * ClrCo=ClrCoSgn="#35978F" 
        (bluegreen)
     * ClrEx=ClrExSgn="#BF812D" 
        (golden)
     * Dash=dashed
     * Ew=prop and 
     * Nw=rel_prop

  * dashed(Dash=dashed)
     defines the edge style for edges with pval >= 0.05 (alt: solid)

  * directed(Dir=false)
     default removes direction from BN

  * edge_metrics(Emetrs=_)
     if present, returns the edge metrics used, in a list of X-Y-po(Pval,Odds)

  * edge_width(Ew=prop)
     defines the edge width for edges (_prop_ or _fix_)

  * node_pen_width(Nw=rel_prop)
     defines the edge width for edges (_rel_prop_ or _fix_)

  * postfix(Psfx='')
    postfix for the output file (postfixing input). Defaults to '' if Theme == org
    and to the Theme otherwise

Requires pack(real) with RcolorBrewer installed.

@author nicos angelopoulos
@version  0.1 2015/11/23
@version  0.2 2018/2/16,    better themes support, and better bic theme
@see fisher.test() in R

*/

gbn_fisher_net( DatF, BnF, DotF ) :-
    gbn_fisher_net( DatF, BnF, DotF, [] ).

gbn_fisher_net( DatF, BnF, DotF, Args ) :-
    debug( gbn(fisher_net), 'Going fishing with: ~p, ~p and ~p', [DatF,BnF,DotF] ),
    options_append( gbn_fisher_net, Args, Opts ),
    options( adjust(Adj), Opts ),
    gbn_fisher_metrics( DatF, Adj, PrvDataPl, IntDfRv, PvalsIntRv, PvalsRv, OddsRv ),

    options( clr_co(Cco), Opts ),
    options( clr_co_signif(Ccs), Opts ),
    options( clr_ex(Cex), Opts ),
    options( clr_ex_signif(Ces), Opts ),
    options( dashed(Dash), Opts ),
    options( edge_width(Ew), Opts ),
    options( node_pen_width(Nw), Opts ),

    gbn_fisher_nodes_attributes( Nw, PrvDataPl, NAttrs ),
    % debug( gbn(fisher_net), 'Nodes attrs: ~w', [NAttrs] ),

    gbn_term( BnF, Bn, _ ),
    debug( gbn(fisher_net), 'Bn read from file: ~p, is: ~w', [BnF,Bn] ),
    maplist( gbn_fisher_family_metrics(PvalsRv,OddsRv), Bn, MeTripsPrv ),
    flatten( MeTripsPrv, MeTrips ),
    memberchk( edge_metrics(MeTrips), Opts ),

    % maplist( gbn_fisher_family_dot_edge_attrs(Style,EwDw,Clr,UpClr,JsDw,JsUp,Dash), Trips, FamsAttrs ),
    maplist( gbn_fisher_family_dot_edge_attrs(Ew,Ces,Ccs,Cex,Cco,Dash), MeTrips, FamsAttrs ), %fixme: double/triple check colours ! 18.02.19
    % EAGoal = gbn_fisher_family_dot_edges(PvalsRv,OddsRv,Ces,Ccs,Cex,Cco,Dash,Ew), 
    % maplist( EAGoal, Bn, EAttrsNest ),
    % flatten( EAttrsNest, EAttrs ),

    gbn_fisher_dot_file( DotF, fclr, BnF ),
    options( postfix(Psfx), Opts ),
    os_postfix( Psfx, DotF, PsfxDotF ),
    os_ext( _, DotS, PsfxDotF ),
    options( bground(ClrBg), Opts ),
    (atomic(ClrBg) -> atom_string(ClrBg,ClrBgStg); ClrBgStg=ClrBg),
    options( directed(Gdir), Opts ),
    ( Gdir == false -> Gtype = graph; Gtype = digraph ),
    Dbns = [ output(svg), dot_file(PsfxDotF), output_stem(DotS), rmv(false),
             graph(bgcolor(ClrBgStg)), node_shape(box), node_style(rounded),
             nodes_dichromatic(false,false),nodes_attrs(NAttrs),edges_attrs(FamsAttrs),
             colour(bnw),type(Gtype)
           ],
    append( Opts,  Dbns, Cohed ),
    disp_bn( Bn, Cohed ),

    os_ext( _, csv, PsfxDotF, OutCsvF ),
    findall( row(X,Y,Pva,Odv), member(X-Y-po(Pva,Odv),MeTrips), CsvRows ),
    mtx( OutCsvF, CsvRows ),

    GbnG = gbn_fisher_signif_dot_edges(Ew,PvalsRv,OddsRv,Ces,Ccs,Bn,Fraph,FAttrsNest),
    os_ext( dat, dot, DatF, FisFPrv ),
    os_postfix( fisher, FisFPrv, FisF ),
    gbn_fisher_if( FisF, GbnG, Fraph, NAttrs, FAttrsNest, Opts ), 
    maplist( <<-, [IntDfRv,PvalsIntRv,PvalsRv,OddsRv] ),
    debug( gbn(fisher_net), 'Done fishing for: ~p', BnF ).

/** gbn_fisher_metrics( +DatF, +BnF, +Adj, -PrvDataPl, -RGbnDf, -RGbnIn, -PvalsRv, -OddsRv ).

*/
gbn_fisher_metrics( DatF, Adj, PrvDataPl, RGbnDf, RGbnIn, PvalsRv, OddsRv ) :-
    csv_read_file( DatF, PrvDataPl, [separator(0' )] ),
    length( PrvDataPl, DataNofRows ),
    PrvDataPl = [DataPlHdr,_|TDataPl],
    functor( DataPlHdr, _, Arity ),
    debug( gbn(fisher_net), 'Data read from: ~p, has ~d rows and ~d columns', [DatF,DataNofRows,Arity] ),
    DataPl = [DataPlHdr|TDataPl], % fixme: debugging only???
    assert( data_fishing(DataPl) ),
    RGbnDf = gbn_df,
    RGbnIn = gbn_inters,
    OddsRv = gbn_odds,
    gbn_mtx_df( DataPl, RGbnDf ),
    % cross_column_fisher_test( gbnf, gbnf_inters, gbnf_odds, Lods ),
    NoClm  <- ncol(RGbnDf),
    RGbnIn <- matrix( ncol=NoClm, nrow=NoClm ),
    OddsRv <- matrix( ncol=NoClm, nrow=NoClm ),
    cross_column_fisher_test( 1, 1, NoClm, RGbnDf, RGbnIn, OddsRv ),
    rownames( RGbnIn ) <- colnames( RGbnDf ),
    colnames( RGbnIn ) <- colnames( RGbnDf ),
    os_ext( bh, RGbnIn, PvalsRv ),
    % dot_ext( RGbnIn, bh, RGbnBh ),
    % Intersbh <- interactions,
    PvalsRv <- RGbnIn,
    ( Adj == false ->
        true
        ;
        % fixme: change the corrected/adjusted R variable name to reflect that non BH corrections might have been applied
        % 
        % atomic_list_concat( ['p.adjust(gbn_inters.bh[lower.tri(gbn_inters.bh)], method="',Adj,'")'], '', Radjust ),
        % 'gbn_inters.bh[lower.tri(gbn_inters.bh)]' <- Radjust,
        atom_string( Adj, AdjSt ),
        'gbn_inters.bh[lower.tri(gbn_inters.bh)]' <- p.adjust('gbn_inters.bh'['lower.tri'('gbn_inters.bh')], method=+AdjSt),
        % RGbnBh[lower.tri(RGbnBh)] <- p.adjust(10^-abs(RGbnBh[lower.tri(RGbnBh)]), method="BH"),
        % 'gbn_inters.bh[upper.tri(gbn_inters.bh)]' <- 'p.adjust(gbn_inters.bh[upper.tri(gbn_inters.bh)], method="BH")'
        'gbn_inters.bh[upper.tri(gbn_inters.bh)]' <- p.adjust('gbn_inters.bh'['upper.tri'('gbn_inters.bh')], method=AdjSt)
    ).

gbn_fisher_family_metrics( BHsR, Odds, X-Ys, Trips ) :-
    RwNms <- rownames(BHsR),
    once( nth1(Xp,RwNms,X) ),
    findall( X-Y-po(Pval,Oval), (
                                    member(Y,Ys), nth1(Yp,RwNms,Y),
                                    Pval <- BHsR[Xp,Yp],
                                    Oval <- Odds[Xp,Yp]
                                ),
                                    Trips ).

gbn_fisher_family_dot_edge_attrs( Ew, DwClr, UpClr, JsDw, JsUp, Dash, X-Y-po(Pval,Oval), FamAttrs ) :-
    gbn_fisher_edge_width( Ew, Oval, Pw ),
    ( Pval < 0.05 ->
        Style = solid,
        ( Oval > 1 -> Clr = UpClr ; Clr = DwClr )
        ;
        Style = Dash,
        ( Oval > 1 -> Clr = JsUp ; Clr = JsDw  )
    ),
    FamAttrs = X-Y-[color(Clr),style(Style),penwidth(Pw)].

gbn_fisher_signif_dot_edges( Ew, BHsR, Odds, DwClr, UpClr, Bn, Graph, NdAttrs ) :-
    RwNms <- rownames(BHsR),
    G = gbn_fisher_node_dot_edges( Ew, BHsR, Odds, DwClr, UpClr, RwNms ),
    maplist( G, Bn, Graph, NdAttrs ).

gbn_fisher_node_dot_edges( Ew, BHsR, Odds, DwClr, UpClr, RwNms, Nd-_, Nd-Fam, Attrs ) :-
    once( nth1(Xp,RwNms,Nd) ),
    findall( Nd-Y-[color(Clr),penwidth(Pw)], (  nth1(Yp,RwNms,Y), Yp =\= Xp,
                            Pval <- BHsR[Xp,Yp],
                            Oval <- Odds[Xp,Yp],
                            Pval < 0.05,
                            gbn_fisher_edge_width(Ew,Oval,Pw),
                            ( Oval > 1 -> 
                                Clr = UpClr
                                ;
                                Clr = DwClr
                            )
                          ),
                            Attrs ),
    findall( Y, member(_-Y-_,Attrs), FamL ),
    sort( FamL, Fam ).

gbn_fisher_if( FisF, _GbnG, _Fraph, _NAttrs, _FAttrsNest, _Opts ) :-
    exists_file( FisF ),
    !,
    debug( gbn(fisher_net), 'Skipping existing Fisher file: ~p', FisF ).
gbn_fisher_if( FisF, Goal, Fraph, NAttrs, FAttrsNest, Opts ) :-
    % os_ext( _, FisS, PsfxFisF ),
    call( Goal ),
    flatten( FAttrsNest, FAttrs ),
    Fbns = [ output(svg), dot_file(FisF), rmv(false),
            graph(bgcolor("#FFFFFF")), node_shape(box), node_style(rounded),
            nodes_dichromatic(false,false),nodes_attrs(NAttrs),edges_attrs(FAttrs),colour(bnw)
           ],
    append( Opts,  Fbns, Fohed ),
    disp_bn( Fraph, Fohed ).

gbn_fisher_dot_file( DotF, Psfix, BnF ) :-
    var( DotF ),
    !,
    os_ext( _, dot, BnF, IntDotF ),
    os_postfix( Psfix, IntDotF, DotF ).
gbn_fisher_dot_file( _DotF, _Psfix, _BnF ).

gbn_fisher_nodes_attributes( fix, _PrvDataPl, [] ).
gbn_fisher_nodes_attributes( rel_prop, Mtx, NAttrs ) :-
    mtx_value_column_frequencies( Mtx, 1, Freqs ),
    kv_decompose( Freqs, Lbls, Times ),
    list_proportions( Times, Propos, to_range(r(1,4)) ),
    findall( Lbl-[penwidth(Nw)], di_member(Lbls,Propos,Lbl,Nw), NAttrs ).

% fixme: goes to stoics_lib ???
di_member( [H1|_T1], [H2|_T2], H1, H2 ).
di_member( [_H1|T1], [_H2|T2], E1, E2 ) :-
    di_member( T1, T2, E1, E2 ).

% fixme: add edges that are significant but not in Ys ? 
% fixme: add parameters to prop,... controlling overall width and factors ???
gbn_fisher_edge_width( prop, Oval, Pw ) :-
    ( number(Oval) -> atom_number(OvalAtm,Oval); OvalAtm = Oval ), % should go straight to if part below if not a number
    % fixme: is the ONan test too general ?
    ( (atom_concat(_,'NaN',OvalAtm);Oval=='1.0Inf';Oval=:=0;Oval=:= 1.0Inf ) -> 
        Pw is 4
        ;
        ( Oval > 1 ->
            Pw is 1.8 + min(3,Oval/4)
            ;
            % Pw is 1 + min(3,log(min(2*(1/Oval),10)))
            % golden: mutual.exclusivity
            Pw is 1.8 + min(3,(1/Oval)/4)
            % Pw is max( 1, 4 - (6 * Oval) )
        )
    ).
gbn_fisher_edge_width( fix, _Oval, 1 ).

/*
gbn_fisher_gated_edges( [], Fraph, Fraph, [] ).
gbn_fisher_gated_edges( [Fam,Stats|T], Fraph, GtGraph, GtAttrs ) :-
    Fam =.. [row,Ch|Pas],
    arg( 6, Stats, Gate ),
    gbn_gate_simplify( Gate, Sate ),
    % ( select(Ch-Pas,Fraph,RFraph) -> true; throw(missing_family(Ch,Pas)) ),
    gbn_gate_family_edges( Pas, Ch, Pas, Gate, GtGraph, TGtGraph, GtAttrs, TGtAttrs ),
    gbn_fisher_gated_edges( T, RFraph, TGtGraph, TGtAttrs ).

gbn_gate_family_edges( Ch, Pas, Gate, GtGraph, TGtGraph, GtAttrs, TGtAttrs ) :-
*/

gbn_gate_simplify( o(A,o(B,C)), o(A,B,C) ) :- !.
gbn_gate_simplify( a(A,a(B,C)), a(A,B,C) ) :- !.
gbn_gate_simplify( x(A,x(B,C)), x(A,B,C) ) :- !.
gbn_gate_simplify( A, A ).
