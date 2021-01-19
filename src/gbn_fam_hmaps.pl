
:- lib(multi_cow_plot/2).
:- lib(mtx_mut_hmap/2).
    % were in pack(sanger), now local

% :- lib(mtx).
% :- lib(real).
% :- lib(options).

% :- lib( stoics_lib:kv_decompose/3 ). % loaded at the top

gbn_fam_hmaps_defaults( Defs ) :-
    Defs = [ cellwidth(1),
             cellheight(12),
             col_mut("#CB181D"),
             col_wt("#08519C"),
             outputs(png),
             x11(X11B)
           ],
    ( getenv('SSH_TTY',_Shy) -> X11B = false; X11B = true ).

/** gbn_fam_hmaps( +Opts ).

Create family heatmaps for all families and all Bns in directory Dir (default is the current directory).

Every file ending in .bn is taken to a (Gobnilp generated) BN. There should be a single .dat file in the directory,<br>
that is taken to hold the data used to create all BNs (Gobnilp default format).

Opts:
  * cellwidth(Cw=1)
     cell width
  * cellheight(Ch=12)
     cell height
  * col_mut(ClrM="#CB181D")
      colour of mutation
  * col_wt(ClrW="#08519C")
      colour of wild type
  * dir(Dir='.')
      working directory
  * outputs(Outs=png)
      output format(s)
  * x11(X11B)
     defaults to false if SSH_TTY is defined OS variable, and true otherwise; passed to mtx_mut_hmap/2

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
	% Width = 960,
    % Dlists = [Fdl|_],
    % length( Fdl, Ldl ),
    % options( cellwidth(Cwidth), Opts ),
    % options( cellheight(Cheight), Opts ),
    % Width is  ( ( ( Ldl // 200 ) + 1 ) * 200 *Cwidth ) + 200,  % fixme: Cwidth
	% Width = 1200,  % use for mpn ( make this dependendant on number of columns !?
    os_path( Dir, GoBnF, GoBnP ),
    debug( gbn(fam_hmaps), 'Family heatmaps for file: ~p', GoBnP ),
	gbn_term( GoBnP, GoBnPrv, _Ado ),
	% findall( BnNode, member(Nd-_,GoBn), BnNodes ),
	os_ext( bn, Stem, GoBnP ),
	atomic_list_concat( [Stem,fams], '_', FamsD ),
	atomic_list_concat( [Stem,prns], '_', PrnsD ),
	os_make_path( FamsD ),
	os_make_path( PrnsD ),
    options( outputs(Outs), Opts ),
    gbn_fam_hmaps_order_net( GoBnPrv, GoBn ),
	findall( [LeadRow,BestRow]-(MmhN-PltN), (
                  nth1( N, GoBn, Node-Pas ),
                  % debug( gbn(fam_hmaps), 'Node: ~w, parents: ~w', [Node,Pas] ),
				  % member( Node-Pas, GoBn ),
				  findall( Ch, (member(Ch-ChPas,GoBn),memberchk(Node,ChPas)), Chs ),
				  \+ (Pas == [], Chs == [] ),
				  append( Pas, [Node|Chs], Family ),
                  % length( Family, FaLen ),
                  % FamHeight is 200 + ( (FaLen - 1) * 20 ),
				  findall( Row, ( member(ANode,Family), memberchk([ANode|Clm],Dlists),
				                  Row =.. [row,ANode|Clm]
				                ),
							 	AllRows ),
				  % os_ext(svg,ANode,RelSvgF ),
                  %
                  % parents only: 
				  os_path(FamsD,Node,NodeStem),
                  % os_postfix(gates,NodeStem,GatesF,ext(csv)),
                  % test taht Dlists is accepted in below call
                  os_ext( csv, NodeStem, NodesCsvF ),
                  % 
                  % PostG = <-( legend(0,0,c("a","b"),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red")) ),
                  % take 3 !
                  % HmOpts = [plot(Plot),rvar(MmhN),stem(NodeStem),outputs(png)],    % outputs(png(width=Width,height=FamHeight)),
                  % FaOpts = [plot(Plot),stem(NodeStem),outputs(png)],
                  options( x11(X11), Opts ),
                  FaOpts = [x11(X11),stem(NodeStem),outputs(Outs)],
                  ( current_predicate(user:display_var_as/2) ->
                        maplist( change_row_name, AllRows, TmpRows )
                        ;
                        AllRows = TmpRows
                  ),
                  mtx_mut_hmap( TmpRows, FaOpts ),
				  Pas \== [],
                  % fixme: Pas is probably already sorted...
                  sort( Pas, OPas ),
                  atomic_list_concat( [Node|OPas], '.', PasNodeBase ),
				  os_path(PrnsD,PasNodeBase,PasNodeStem),
				  append( Pas, [Node], PasNodeL ),
				  findall( Row, ( member(ANode,PasNodeL), memberchk([ANode|Clm],Dlists),
				                  Row =.. [row,ANode|Clm]
				                ),
							 	PasRows ),
                  atomic_list_concat( [mmh,N], '_', MmhN ),
                  ( current_predicate(user:display_var_as/2) ->
                        maplist( change_row_name, PasRows, TmpPasRows )
                        ;
                        PasRows = TmpPasRows
                  ),
                  ( memberchk(x11(X11),Opts) ->
                    PaOpts = [x11(X11),plot(Plot),rvar(MmhN),stem(PasNodeStem),outputs(Outs)]    % outputs(png(width=Width,height=FamHeight)),
                    ;
                    PaOpts = [plot(Plot),rvar(MmhN),stem(PasNodeStem),outputs(Outs)]    % outputs(png(width=Width,height=FamHeight)),
                  ),
                  mtx_mut_hmap( TmpPasRows, PaOpts ),
                  atomic_list_concat( [plt,N], '_', PltN ),
                  PltN <- Plot,
                  gbn_family_gates(Node,Pas,Dlists,Gatrix,Opts),
                  %
                  ( Pas = [_] ->
                    LeadRow = [], BestRow = []
                    ;
                    os_postfix( gates, NodesCsvF, GatesF ),
                    mtx( GatesF, Gatrix ),
                    Gatrix = [BestRow|_], 
                    LeadRow=.. [row,Node|Pas]
                  )
				), 
					NestPrs ),
    kv_decompose( NestPrs, Nest, Nrs ),
    kv_decompose( Nrs, Mtvs, Plvs ),
    flatten( Nest, LeadsBests ),
    % os_ext( csv, PrnsD, PrnsCsvF ),

    atomic_list_concat( [Stem,gates,best], '_', GatesStemBF ),
    os_ext( csv, GatesStemBF, GatesBF ),
    % os_postfix( gates_best, PrnsCsvF, GatesBF ),
    % os_postfix( multi_prns, PrnsCsvF, MultiFamsF ),
    % os_postfix( multi_prns, PrnsCsvF, MultiFamsF ),
	atomic_list_concat( [Stem,multi,prns], '_', MultiPrnsStem ),  % was MultiFamStem
    % Stem
    % os_ext( _, MultiFamStem, MultiFamsF ),
    ( Plvs == [] ->
        debug( gbn(fam_hmaps), 'empty list for multi-plot in gbn_fam_hmaps/1', true )
        ;
        % fixme: this will breake if Outs is not an atomic extension...
        % write( plvs(Plvs) ), nl,
        multi_cow_plot( Plvs, [stem(MultiPrnsStem),ext(Outs),labels(lower)|Opts] )
    ),
    maplist( r_remove, Mtvs ),
    maplist( r_remove, Plvs ),
    csv_write_file( GatesBF, LeadsBests, [match_arity(false)] ).

gbn_fam_hmaps_order_net( GoBn, Order ) :-
    current_predicate( gbn:hmap_family_order/1 ),
    !,
    findall( Ch-Pa, (gbn:hmap_family_order(Ch),memberchk(Ch-Pa,GoBn)), Left ),
    findall( Ch-Pa, (member(Ch-Pa,GoBn),\+ member(Ch-Pa,Left)), Right ),
    append( Left, Right, Order ).
    % write( append(Left,Right,Order) ), nl.
gbn_fam_hmaps_order_net( GoBn, GoBn ).

change_row_name( RowIn, RowOut ) :-
    RowIn =.. [Tname,Rname|Args],
    ( user:display_var_as(Rname,Nname) -> true; Rname = Nname ),
    RowOut =.. [Tname,Nname|Args].
