
:- use_module(library(debug)).
:- use_module(library(lib)).

:- lib(gbn).
:- lib(real).
:- lib(b_real).
:- lib(os_lib).
:- lib(debug_call).

:- debug(aml).

/** aml_edge_penalty.

Build an AML edge penalty plot, where for selected number of variables min60,
and varying ε (edge_penalty) variables we plot the number of 
edges on the resulting BN.

@author nicos angelopoulos
@version  0:1 2020/07/25

*/
aml_edge_penalty :-
        debug_call( aml, start, var_sel_plot ),
    % aml_edge_penalty_run,
    aml_ep_plot,
        debug_call( aml, end, var_sel_plot ).

aml_edge_penalty_run :-
        debug_call( aml, start, var_sel_run ),
    DatF = 'data/aml/aml_min60.dat',
    OsOdir = aml_min60_penalty,
        debug_call( aml, start, file(DatF) ),
    numlist( 1, 20, Es ),
    GBNOpts = [ data(DatF),
                multiple(edge_penalty,Es),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debug_call( aml, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debug_call( aml, start, fam_hmaps ),
    gbn_fam_hmaps( dir(OsOdir) ),   % only for binaries
        debug_call( aml, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debug_call( aml, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debug_call( aml, end, var_sel_run ).

aml_ep_plot :-
    Dir = aml_min60_penalty,
    os_sel( os_files, [ext(dot),postfix(fclr)], Files, dir(Dir) ),
    maplist( aml_file_edge_count_kv(Dir), Files, KVs ),
    sort( KVs, OrdKVs ),
    maplist( writeln, OrdKVs ),
    % FClrs = ["darkolivegreen"],
    FClrs = ["#779ECB"],
    Opts = [ flip(false), 
             geom_bar(empty),
             labels('Edge penalty value','Number of Edges','Effect of edge penalty on network size.'),
             fill_colours(FClrs),
             panel_theme(axes)
           ],
    gg_bar_plot( OrdKVs, Opts ),
    <- ggsave( "pdfs/aml_ep_plot.pdf", width=14 ).

aml_file_edge_count_kv( Dir, File, E-Count ) :-
    directory_file_path( Dir, File, Path ),
    edge_count( Path, Count ),
    atomic_list_concat( [_,Psfx], '-e', File ),
    atomic_list_concat( [Etm|_], '_', Psfx ),
    atom_number( Etm, E ).

/** edge_count( +DotF, -EdgesCount ).

Count number of edges in a Dot File.

==
?- edge_count( 'aml_min60-e10_fclr.dot', Count ).
Count = 24.
==

@author nicos angelopoulos
@version  0:1 2020/07/26

*/
edge_count( DotF, Edges ) :-
    gbn:dot_read( DotF, dot(_,_,_,Terms) ),
    findall( 1, member(edge(_,__),Terms), Ones ),
    sumlist( Ones, Edges ).
