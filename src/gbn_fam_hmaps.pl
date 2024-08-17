
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
             multi_prns(true),
             mut_hmap_iface(pl),
             outputs(png),
             plot_fams(true),
             plot_prns(true),
             x11(X11B)
           ],
    ( getenv('SSH_TTY',_Shy) -> X11B = false; X11B = true ).

/** gbn_fam_hmaps( +Opts ).

Create family and parental heatmaps for all families and all Bns in directory Dir (default is the current directory).

Every file ending in .bn is taken to be a (Gobnilp generated) BN. 
There should be a single .dat file in the directory, which is taken to hold
the data used to create all BNs (Gobnilp default format).

The predicate creates 2 types of heatmaps: parental and family based. Family ones includes the children of a node 
in addition to its parents. Currently the parental heatmaps are also summarised in a single multiplot.

Opts:
  * as_mutational(Bin=true)
    whether the data contains mutations<br> (when false arbitary integer values are assumed)
  * cellwidth(Cw=1)
    cell width
  * cellheight(Ch=12)
    cell height
  * col_hmap(ClrH=[])
    list of value-"Colour" for each value in mtx (can contain non mtx values).<br>
    If this is given it overrides ClrM and ClrW.
  * col_mut(ClrM="#CB181D")
    colour of mutation
  * col_wt(ClrW="#08519C")
    colour of wild type
  * debug(Dbg=false)
    use _gbn(fam_hmaps)_ for high level messages and _gbn(fam_hmaps_fine)_ for finner grain.<br>
    The latter also sets the former. Only the first debug() is observed.
  * dir(Dir='.')
    working directory
  * mut_hmap_iface(Ifc=pl)
    Interface for calling mtx_mut_hmap/1,2. Use _os_ for an external call to upsh via shell/1.
  * multi_prns(MuPrns=true)
    whether to produce a multi parental plot
  * outputs(Outs=png)
    output format(s) - propagates to multi_cow_plot/2 as ext()
  * plot_fams(PtFams=true)
    enables family heatmap plotting
  * plot_prns(PtFams=true)
    enables parents heatmap plotting
  * x11(X11B)
    defaults to false if SSH_TTY is defined OS variable, and true otherwise.<br>
    Passed to mtx_mut_hmap/2 and multi_cow_plot/2

Options are also passed to gbn_family_gates/5 and multi_cow_plot/2.

@author nicos angelopoulos
@version  0.2 2018/02/21
@tbd add a simple example in examples/
@tbd multi family plot ?
@tbd add token to allow multi outputs in same directory

*/
gbn_fam_hmaps( Args ) :-
    options_append( gbn_fam_hmaps, Args, Opts ),
    options( dir(Dir), Opts ),
    debug_chain( fam_hmaps_fine, fam_hmaps ),
    debuc( gbn(fam_hmaps), 'Starting family and parental heatmaps in dir: ~p', [Dir] ),
    debuc( gbn(fam_hmaps_fine), 'Starting (finely) family and parental heatmaps in dir: ~p', [Dir] ),
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
    debuc( gbn(fam_hmaps), 'Family heatmaps for file: ~p', [GoBnP] ),
    gbn_term( GoBnP, GoBnPrv, _Ado ),
    os_ext( bn, Stem, GoBnP ),
    atomic_list_concat( [Stem,fams], '_', FamsD ),
    atomic_list_concat( [Stem,prns], '_', PrnsD ),
    os_make_path( FamsD ),
    os_make_path( PrnsD ),
    options( outputs(Outs), Opts ),
    gbn_fam_hmaps_order_net( GoBnPrv, GoBn ),
    options( plot_fams(PtFams), Opts ),
    options( plot_prns(PtPrns), Opts ),
    gbn_fam_hmaps_plots( GoBn, 1, GoBn, Dlists, PtFams/PtPrns, FamsD, PrnsD, Nest, Nrs, Opts ),
    % kv_decompose( NestPrs, Nest, Nrs ),
    kv_decompose( Nrs, Mtvs, Plvs ),
    flatten( Nest, LeadsBests ),
    atomic_list_concat( [Stem,gates,best], '_', GatesStemBF ),
    os_ext( csv, GatesStemBF, GatesBF ),
    options( multi_prns(MuPrns), Opts ),
    ( Plvs == [] ->
        debuc( gbn(fam_hmaps), 'empty list for multi-plot in gbn_fam_hmaps/1', [] )
        ;
        ( MuPrns == false ->
               debuc( gbn(fam_hmaps), 'Skipping creation of multi parental plot due to flag value.', [] )
               ;
               atomic_list_concat( [Stem,multi,prns], '_', MultiPrnsStem ),
               % fixme: this will break if Outs is not an atomic extension...
               multi_cow_plot( Plvs, [stem(MultiPrnsStem),ext(Outs),labels(lower)|Opts] )
        )
    ),
    maplist( r_remove, Mtvs ),
    maplist( r_remove, Plvs ),
    csv_write_file( GatesBF, LeadsBests, [match_arity(false)] ).

gbn_fam_hmaps_plots( [], _I, _GoBn, _Dlists, _Pts, _FamsD, _PrnsD, [], [], _Opts ).
gbn_fam_hmaps_plots( [Node-Pas|Bn], N, GoBn, Dlists, PtFs/PtPs, FamsD, PrnsD, [LeadRow,BestRow|LBRs], MPs, Opts ) :-
     debuc( gbn(fam_hmaps_fine), 'Node: ~w, with parents: ~w', [Node,Pas] ),
     findall( Ch, (member(Ch-ChPas,GoBn),memberchk(Node,ChPas)), Chs ),
     Pas \== [],
     \+ (Pas == [], Chs == [] ),
     !,
     options( mut_hmap_iface(Ifc), Opts ),
     ( PtFs == true ->
          flatten( [Pas,Chs,Node], Family ),  % 24.08.16: makes id-ing the central node easier [was append(Pas,[Node|Chs],Family)]
                                              % fixme: the order is changed later...
          debuc( gbn(fam_hmaps_fine), 'Getting rows for family: ~w', [Family] ),
          gbn_fam_hmaps_rows( Dlists, Family, TmpRows ),
          debuc( gbn(fam_hmaps_fine), length, fam_rows/TmpRows ),
          % findall(Row1, (member(ANode1,Family),memberchk([ANode1|Clm1],Dlists),Row1 =.. [row,ANode1|Clm1]), AllRows),
          % parents only: 
          os_path( FamsD, Node, NodeStem ),
          os_ext( csv, NodeStem, NodesCsvF ),
          options( [outputs(Outs),x11(X11)], Opts ),
          options( [as_mutational(Bin),col_wt(ClrW),col_mut(ClrM),col_hmap(ClrH)], Opts ),
          ComOpts = [as_mutational(Bin),col_wt(ClrW),col_mut(ClrM),col_hmap(ClrH)],
          FaOpts = [x11(X11),stem(NodeStem),outputs(Outs)|ComOpts],
          gbn_fam_mut_hmap( Ifc, TmpRows, NodeStem, FaOpts ),
          % mtx_mut_hmap( TmpRows, FaOpts ),
          debuc( gbn(fam_hmaps_fine), 'Finished family heatmap.', [] )
          % fixme: Pas is probably already sorted...
          ;
          debuc( gbn(fam_hmaps_fine), 'Skipping plot for family of: ~w.', [Node] )
     ),
     ( PtPs == true ->
          sort( Pas, OPas ),
          atomic_list_concat( [Node|OPas], '.', PasNodeBase ),
          os_path(PrnsD,PasNodeBase,PasNodeStem),
          append( Pas, [Node], PasNodeL ),
          debuc( gbn(fam_hmaps_fine), 'Getting rows for parents and child: ~w', [PasNodeL] ),
          gbn_fam_hmaps_rows( Dlists, PasNodeL, TmpPasRows ),
          % findall( Row2, (member(ANode2,PasNodeL), memberchk([ANode2|Clm2],Dlists),Row2 =.. [row,ANode2|Clm2]), PasRows ),
          atomic_list_concat( [mmh,N], '_', MmhN ),
          ( memberchk(x11(X11),Opts) ->
               PaOpts = [x11(X11),plot(Plot),rvar(MmhN),stem(PasNodeStem),outputs(Outs)|ComOpts]    % outputs(png(width=Width,height=FamHeight)),
               ;
               PaOpts = [plot(Plot),rvar(MmhN),stem(PasNodeStem),outputs(Outs)|ComOpts]    % outputs(png(width=Width,height=FamHeight)),
          ),
          gbn_fam_mut_hmap( Ifc, TmpPasRows, PasNodeStem, PaOpts ),
          % mtx_mut_hmap( TmpPasRows, PaOpts ),
          atomic_list_concat( [plt,N], '_', PltN ),
          PltN <- Plot,
          MPs = [MmhN-PltN|TMPs]
          ;
          TMPs = MPs,
          debuc( gbn(fam_hmaps_fine), 'Skipping plot for parents of: ~w.', [Node] )
     ),
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
     !,
     garbage_collect,
     gbn_fam_hmaps_plots( Bn, N1, GoBn, Dlists, PtFs/PtPs, FamsD, PrnsD, LBRs, TMPs, Opts ).
% fixme: double check this: (also steadfast Orphan ?, also check if J = I is fine ?
gbn_fam_hmaps_plots( [_Orphan|Bn], I, GoBn, Dlists, Pts, FamsD, PrnsD, LBRs, MPs, Opts ) :-
     J is I + 1,
     gbn_fam_hmaps_plots( Bn, J, GoBn, Dlists, Pts, FamsD, PrnsD, LBRs, MPs, Opts ).

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

gbn_fam_hmaps_rows( [], Rels, Rows ) :-
    ( Rels == [] ->
               Rows = []
               ;
               throw( could_not_find_heatmap_row_for_relatives(Rels) )
    ).
gbn_fam_hmaps_rows( [L|Ls], Rels, Rows ) :-
     L = [ANode|Clm],
     ( select(ANode,Rels,Rems) ->
               ( user:display_var_as(ANode,Dname) -> true; ANode = Dname ),
               Row =.. [row,Dname|Clm],
               Rows = [Row|Tows],
               ( Rems = [] ->
                         Ms = []
                         ;
                         Ls = Ms
               )
               ;
               Ms = Ls, Tows = Rows, Rems = Rels
     ),
     gbn_fam_hmaps_rows( Ms, Rems, Tows ).

gbn_fam_mut_hmap( pl, Rows, _OStem, Opts ) :-
          mtx_mut_hmap( Rows, Opts ),
          !.
gbn_fam_mut_hmap( os, Rows, OStem, Opts ) :-
     os_ext( csv, OStem, CsvStem ),
     os_postfix( data_rows, CsvStem, CsvF ),
     mtx( CsvF, Rows ),
     findall( Atom, (  member(Opt,[mtx(CsvF)|Opts]),
                       (Opt =.. [F,A] -> true; throw(non_uninary_opt_in(Opt))),
                       ( is_list(A) ->
                              gbn_fam_mut_hmap_lst_os( A, AAtm )
                              ;
                              ( string(A) ->
                                        atomic_list_concat( [stg,A], '=', AAtm )
                                        ;
                                        A = AAtm
                              )
                       ),
                       atomic_list_concat( [F,AAtm], '=', Atom )
                    ),
                         Atoms ),
     atomic_list_concat( [upsh,'gbn:mtx_mut_hmap'|Atoms], ' ', Shell ),
     write( shell(Shell) ), nl,
     debuc( gbn(fam_hmaps_fine), 'Shellling: ~w', [Shell] ),
     shell( Shell ).

gbn_fam_mut_hmap_lst_os( [], 'lst=' ).
gbn_fam_mut_hmap_lst_os( [H|T], Atm ) :-
     atomic_list_concat( [H|T], ',', List ),
     atomic_list_concat( [lst,List], '=', Atm ).
