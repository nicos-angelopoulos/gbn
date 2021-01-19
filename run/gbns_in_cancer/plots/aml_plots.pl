
:- use_module(library(debug)).
:- use_module(library(lib)).

:- lib(gbn).
:- lib(mtx).
:- lib(mlu).
:- lib(real).
:- lib(bio_db).
:- lib(stoics_lib).
:- lib(debug_call).

:- [aml_edge_penalty].

:- debug(aml).
:- debug(real).

aml_colour(darkpastelred, '#C23B22').
aml_colour(airforceblue, '#5D8AA8').

/** aml_var_sel.

Build the AML plots. 

1. Population, where for selected ε (edge_penalty) 
and varying number of selected variables we plot the number of 
edges on the resulting BN.

2. Minimum occurances cut-off vs number of variables.

3. Monochrome population plot of 
   (no label or tick for min.mut used line) 

4. Occurance coloured version of the above (w/ all labels present and correct)


@author nicos angelopoulos
@version  0:1 2020/07/25

*/
aml_plots :-
        debug_call( aml, start, var_sel_plot ),
    Cuts = ['00','05','10','20','30','40','50','60','70','80'],
    % fixme: automate so that only run this if not there
    % maplist( aml_run_cut, Cuts ),
    maplist( aml_vs_edge_count_kvs, Cuts, KVs ),
    sort( KVs, OrdKVs ),
    write( ord_kvs(OrdKVs) ), nl,
    FClrs = ["#779ECB"],
    Opts = [ flip(false), 
             geom_bar(empty),
             labels('Cut-off for minimum occurances','Number of Edges','Effect of cut-off on network size.'),
             fill_colours(FClrs),
             panel_theme(axes)
           ],
    gg_bar_plot( OrdKVs, Opts ),
    shell( 'mkdir -p pdfs' ),
    <- ggsave( "pdfs/aml_cut_sel_edge_counts.pdf" ),

    maplist( aml_vs_pop_kvs, Cuts, Pops ),
    write( pops(Pops) ), nl,
    FClrs = ["#779ECB"],
    Ppts = [ flip(false), 
             geom_bar(empty),
             labels('Cut-off for minimum occurances','Number of Variables','Effect of cut-off on number of variables.'),
             fill_colours(FClrs),
             panel_theme(axes)
           ],
    gg_bar_plot( Pops, Ppts ),
    <- ggsave( "pdfs/aml_pops_edge_counts.pdf" ),
    maplist( aml_freqs, ['00'], _ ),
        debug_call( aml, end, var_sel_plot ).

aml_freqs( Cut, _ ) :-
    atomic_list_concat( ['aml_min',Cut], OsOdir ),
    os_ext( dat, OsOdir, DatF ),
    os_path( OsOdir, DatF, DatP ),
    mtx( DatP, MtxIn, sep(0' ) ),
    MtxIn = [Hdr,_ValTot|Rows],
    Mtx = [Hdr|Rows],
    mtx_value_column_frequencies( Mtx, 1, Freqs ),
    write( freqs(Freqs) ), nl,
    findall( V-K, member(K-V,Freqs), Treqs ),
    sort( Treqs, Oreqs ),
    reverse( Oreqs, Rreqs ),
    findall( V-K, member(K-V,Rreqs), Preqs ),
    os_ext( pdf, OsOdir, PdfF ),
    os_path( pdfs, PdfF, PdfP ),
    % mlu_frequency_plot( Freqs, [interface(gg_bar),output(pdf(+PdfP))] ).
    atomic_list_concat( ['Number of patients per driver event for AML.'], '', Main ),
    Mlu = [ interface(barplot), outputs([pdf(file= +PdfP, width=14, height=7)]),
            main= +Main, 
            xlab= "Drivers",
            ylab= "Patients",
            las=2, cex.names=0.5,
            post_call( <-( abline(h=60,col="#C23B22",lwd=2,lty=2)) )
          ],
    mlu_frequency_plot( Preqs, Mlu ),

    % aml_min00_muts_by_pnms.pdf :
    aml_add_mutations_column( Mtx, Stx ),
    aml_colour( airforceblue, AFBclr ),
    VLine = geom_segment(aes(x= -0.01, y = 60, xend = 84, yend = 60), color= +AFBclr, linetype="dashed", size=0.4),
    % Vtext = annotate("text", x= -3, y= 60, hjust=0, vjust=0, label= "60"),
    % tmp <- 'data.frame'(tmpx=c(1,1),tmpy=c(2,3)),
    YTicks = scale_y_continuous(breaks = c(60,100,200,300,400)),
    YTheme = theme(axis.text.y = element_text(color = c(+AFBclr, "black", "black", "black", "black"), size=c(11,9,9,9,9) ),
         axis.ticks.y = element_line(color = c(+AFBclr, "black", "black", "black", "black"),
                          size = c(1,1,1,1,1))),
    os_postfix( muts_by_pnms, PdfP, ByPnmsPdfP ),
    os_ext( pdf, ByPnmsPdfP, ByPnmsPdfS ),
    Pnms = [panel_theme(axes),gg_terms([VLine,YTicks,YTheme]),x_axis_colour(x_axis_aml)],
    gbn:gg_muts_by_pnms( Stx, ByPnmsPdfS, Pnms ).

x_axis_aml( KVs, Clrs, Leg ) :-
    % write( here(A,B,C) ), nl,
    kv_decompose( KVs, Ks, _Vs ),
    partition( hgnc_symbol_start, Ks, Symbs, _Nons ),
    Left = "#3B444B", % Left = "Arsenic",
    % Right = "#B2BEB5",  % Ash Grey
    Right = "#87A96B", % Asparagus
    split_colours( Symbs, Ks, Left, Right, Clrs ),
    sort( ['Point_mutations'-Left,'Complex_events'-Right], Ord ), 
    Leg = bottom( 'Event_type', 8, 1, Ord ).

split_colours( [], Ks, _Left, Right, Clrs ) :-
    findall( Right, member(_,Ks), Clrs ).
split_colours( [Symb|Symbs], [K|Ks], Left, Right, [Clr|Clrs] ) :-
    ( Symb == K ->
        Rymbs = Symbs,
        Clr = Left
        ;
        Rymbs = [Symb|Symbs],
        Clr = Right
    ),
    split_colours( Rymbs, Ks, Left, Right, Clrs ).

hgnc_symbol_start( Symb ) :-
    hgnc_symbol( Symb ),
    !.
hgnc_symbol_start( Atom ) :-
    atomic_list_concat( [Symb|_], '_', Atom ),
    hgnc_symbol( Symb ),
    !.

aml_add_mutations_column( [Hdr|Rows], [Ndr|Nows] ) :-
    Hdr =.. [Name|Args],
    append( Args, ['Mutations'], Nrgs ),
    Ndr =.. [Name|Nrgs],
    maplist( aml_add_mutations_column_row, Rows, Nows ).

aml_add_mutations_column_row( Row, Now ) :-
    Row =.. [Name|Muts], 
    sumlist( Muts, Tot ),
    append( Muts, [Tot], Nargs ),
    Now =.. [Name|Nargs].

aml_add_sample_column( [Hdr|Rows], [Ndr|Nows] ) :-
    Hdr =.. [Name|Args],
    Ndr =.. [Name,'Sample'|Args],
    aml_add_sample_row( Rows, 1, Nows ).

aml_add_sample_row( [], _I, [] ).
aml_add_sample_row( [R|Rs], I, [N|Ns] ) :-
    R =.. [Name|Args],
    N =.. [Name,I|Args],
    J is I + 1,
    aml_add_sample_row( Rs, J, Ns ).

aml_vs_pop_kvs( Cut, KV ) :-
    atomic_list_concat( ['aml_min',Cut], OsOdir ),
    % aml_min05/aml_min05.dat
    os_ext( dat, OsOdir, DatF ),
    os_path( OsOdir, DatF, DatP ),
    os_ext( dat, OsOdir, DatF ),
    os_path( OsOdir, DatF, DatP ),
    mtx( DatP, Mtx, sep(0' ) ),
    mtx_dims( Mtx, _NofRows, NofVars ),
    KV = Cut-NofVars.

aml_vs_edge_count_kvs( Cut, KV ) :-
    atomic_list_concat( ['aml_min',Cut], OsOdir ),
    % aml_min05/aml_min05.dot
    os_ext( dot, OsOdir, DotF ),
    os_path( OsOdir, DotF, DotP ),
    edge_count( DotP, Count ),
    KV = Cut-Count.

aml_run_cut( Cut ) :-
        debug_call( aml, start, running_cut(Cut) ),
    atomic_list_concat( ['data/aml/aml_min',Cut,'.dat'], DatF ),
    atomic_list_concat( ['aml_min',Cut], OsOdir ),
        debug_call( aml, start, file(DatF) ),
    E = 7,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
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
        debug_call( aml, end, ran_cut(Cut) ).
