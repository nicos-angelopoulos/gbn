
:- lib(multi_cow_plot/2).
:- lib(mtx_mut_hmap/2).
    % were in pack(sanger), now local

% :- lib(mtx).
% :- lib(real).
% :- lib(options).

% :- lib( stoics_lib:kv_decompose/3 ). % loaded at the top

gbn_fam_hmaps_defaults( Defs ) :-
    Defs = [ as_mutational(true),
             cellwidth(1),
             cellheight(12),
             col_hmap([]),
             col_mut("#CB181D"),
             col_wt("#08519C"),
             outputs(png),
             x11(X11B)
           ],
    ( getenv('SSH_TTY',_Shy) -> X11B = false; X11B = true ).

/** gbn_fam_hmaps( +Opts ).

Create family heatmaps for all families and all Bns in directory Dir (default is the current directory).

Every file ending in .bn is taken to be a (Gobnilp generated) BN. 
There should be a single .dat file in the directory,
which is taken to hold the data used to create all BNs (Gobnilp default format).

Opts:
  * as_mutational(Bin=true)
    whether the data contains mutations (when false arbitary integer values are assumed)
  * cellwidth(Cw=1)
    cell width
  * cellheight(Ch=12)
    cell height
  * col_hmap(ClrH=[])
    list of value-"Colour" for each value in mtx (can contain non mtx values).
    If this is given it overrides ClrM and ClrW.
  * col_mut(ClrM="#CB181D")
    colour of mutation
  * col_wt(ClrW="#08519C")
    colour of wild type
  * dir(Dir='.')
    working directory
  * outputs(Outs=png)
    output format(s) - propagates to multi_cow_plot/2 as ext()
  * x11(X11B)
    defaults to false if SSH_TTY is defined OS variable, and true otherwise; 
    passed to mtx_mut_hmap/2 and multi_cow_plot/2

Options are also passed to gbn_family_gates/5 and multi_cow_plot/2.

@author nicos angelopoulos
@version  0.2 2018/02/21
@tbd add a simple example in examples/

*/
gbn_fam_hmaps( Args ) :-
    options_append( gbn_fam_hmaps, Args, Opts ),
    options( dir(Dir), Opts ),
    os_sel( os_files, ext(bn), GoBn, dir(Dir) ),
    gbn_res_dir_dat_file( Dir, DatF ),
    os_path( Dir, DatF, DatP ),
    csv_read_file( DatP, Dat, [separator(0' )] ),
    Dat = [Hdr,_Cnts|Data], 
    Rat = [Hdr|Data],
    mtx_lists( Rat, Dlists ),
    maplist( gbn_fam_hmaps_dlists(Dir,Dlists,Opts), GoBn ).

gbn_fam_hmaps_dlists( Dir, Dlists, Opts, GoBnF ) :-
    % Width is  ( ( ( Ldl // 200 ) + 1 ) * 200 *Cwidth ) + 200,  % fixme: Cwidth
    os_path( Dir, GoBnF, GoBnP ),
    debuc( gbn(fam_hmaps), 'Family heatmaps for file: ~p', GoBnP ),
    gbn_term( GoBnP, GoBnPrv, _Ado ),
    os_ext( bn, Stem, GoBnP ),
    atomic_list_concat( [Stem,fams], '_', FamsD ),
    atomic_list_concat( [Stem,prns], '_', PrnsD ),
    os_make_path( FamsD ),
    os_make_path( PrnsD ),
    options( outputs(Outs), Opts ),
    gbn_fam_hmaps_order_net( GoBnPrv, GoBn ),
    gbn_fam_hmaps_plots( GoBn, 1, GoBn, Dlists, FamsD, PrnsD, Nest, Nrs, Opts ),
    % kv_decompose( NestPrs, Nest, Nrs ),
    kv_decompose( Nrs, Mtvs, Plvs ),
    flatten( Nest, LeadsBests ),
    atomic_list_concat( [Stem,gates,best], '_', GatesStemBF ),
    os_ext( csv, GatesStemBF, GatesBF ),
    atomic_list_concat( [Stem,multi,prns], '_', MultiPrnsStem ),  % was MultiFamStem
    ( Plvs == [] ->
        debuc( gbn(fam_hmaps), 'empty list for multi-plot in gbn_fam_hmaps/1', true )
        ;
        % fixme: this will break if Outs is not an atomic extension...
        multi_cow_plot( Plvs, [stem(MultiPrnsStem),ext(Outs),labels(lower)|Opts] )
    ),
    maplist( r_remove, Mtvs ),
    maplist( r_remove, Plvs ),
    csv_write_file( GatesBF, LeadsBests, [match_arity(false)] ).

gbn_fam_hmaps_plots( [], _I, _GoBn, _Dlists, _FamsD, _PrnsD, [], [], _Opts ).
gbn_fam_hmaps_plots( [Node-Pas|Bn], N, GoBn, Dlists, FamsD, PrnsD, [LeadRow,BestRow|LBRs], [MmhN-PltN|MPs], Opts ) :-
     findall( Ch, (member(Ch-ChPas,GoBn),memberchk(Node,ChPas)), Chs ),
     Pas \== [],
     \+ (Pas == [], Chs == [] ),
     !,
     append( Pas, [Node|Chs], Family ),
     findall( Row1, (member(ANode1,Family),memberchk([ANode1|Clm1],Dlists),Row1 =.. [row,ANode1|Clm1]), AllRows ),
     % parents only: 
     os_path( FamsD, Node, NodeStem ),
     os_ext( csv, NodeStem, NodesCsvF ),
     options( [outputs(Outs),x11(X11)], Opts ),
     options( [as_mutational(Bin),col_wt(ClrW),col_mut(ClrM),col_hmap(ClrH)], Opts ),
     ComOpts = [as_mutational(Bin),col_wt(ClrW),col_mut(ClrM),col_hmap(ClrH)],
     FaOpts = [x11(X11),stem(NodeStem),outputs(Outs)|ComOpts],
     ( current_predicate(user:display_var_as/2) ->
                        maplist( change_row_name, AllRows, TmpRows )
                        ;
                        AllRows = TmpRows
     ),
     mtx_mut_hmap( TmpRows, FaOpts ),
     % fixme: Pas is probably already sorted...
     sort( Pas, OPas ),
     atomic_list_concat( [Node|OPas], '.', PasNodeBase ),
     os_path(PrnsD,PasNodeBase,PasNodeStem),
     append( Pas, [Node], PasNodeL ),
     findall( Row2, (member(ANode2,PasNodeL), memberchk([ANode2|Clm2],Dlists),Row2 =.. [row,ANode2|Clm2]), PasRows ),
     atomic_list_concat( [mmh,N], '_', MmhN ),
     ( current_predicate(user:display_var_as/2) ->
          maplist( change_row_name, PasRows, TmpPasRows )
          ;
          PasRows = TmpPasRows
     ),
     ( memberchk(x11(X11),Opts) ->
          PaOpts = [x11(X11),plot(Plot),rvar(MmhN),stem(PasNodeStem),outputs(Outs)|ComOpts]    % outputs(png(width=Width,height=FamHeight)),
          ;
          PaOpts = [plot(Plot),rvar(MmhN),stem(PasNodeStem),outputs(Outs)|ComOpts]    % outputs(png(width=Width,height=FamHeight)),
     ),
     mtx_mut_hmap( TmpPasRows, PaOpts ),
     atomic_list_concat( [plt,N], '_', PltN ),
     PltN <- Plot,
     gbn_family_gates( Node, Pas, Dlists, Gatrix, Opts ),
     %
     ( Pas = [_] ->
          LeadRow = [], BestRow = []
          ;
          os_postfix( gates, NodesCsvF, GatesF ),
          mtx( GatesF, Gatrix ),
          Gatrix = [BestRow|_], 
          LeadRow=.. [row,Node|Pas]
     ),
     N1 is N + 1,
     gbn_fam_hmaps_plots( Bn, N1, GoBn, Dlists, FamsD, PrnsD, LBRs, MPs, Opts ).
% fixme: double check this: (also steadfast Orphan ?, also check if J = I is fine ?
gbn_fam_hmaps_plots( [_Orphan|Bn], I, GoBn, Dlists, FamsD, PrnsD, LBRs, MPs, Opts ) :-
     J is I + 1,
     gbn_fam_hmaps_plots( Bn, J, GoBn, Dlists, FamsD, PrnsD, LBRs, MPs, Opts ).

gbn_fam_hmaps_order_net( GoBn, Order ) :-
    current_predicate( gbn:hmap_family_order/1 ),
    !,
    findall( Ch-Pa, (gbn:hmap_family_order(Ch),memberchk(Ch-Pa,GoBn)), Left ),
    findall( Ch-Pa, (member(Ch-Pa,GoBn),\+ member(Ch-Pa,Left)), Right ),
    append( Left, Right, Order ).
gbn_fam_hmaps_order_net( GoBn, GoBn ).

change_row_name( RowIn, RowOut ) :-
    RowIn =.. [Tname,Rname|Args],
    ( user:display_var_as(Rname,Nname) -> true; Rname = Nname ),
    RowOut =.. [Tname,Nname|Args].
