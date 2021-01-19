
:- ensure_loaded(aml_plots).
:- debug(etc).

/** etc_plots.

Population plots for all non-AML datasets.

*/

etc_plots :-
    pop_plots( mye, 20, 'Myeloma' ),
    pop_plots( coa,  5,  'TCGA colon adenocarcinoma' ),
    pop_plots( gbm,  5,  'glioblastoma' ),
    pop_plots( mpn,  5,  'myeloproliferative neoplasms' ),
    true.

pop_plots( Canc, Cut, Clb ) :-
        debug_call( etc, start, Canc ),
    % atomic_list_concat( [Canc,'_min',Cut], '', Stem ),
    % atomic_list_concat( [Canc,'_min',Cut], '', Stem ),
    os_ext( dat, Canc, DatF ),
    os_path( Canc, DatF, SubP ),
    os_path( data, SubP, DatP ),
    mtx( DatP, MtxIn, sep(0' ) ),
    MtxIn = [Hdr,_ValTot|Rows],
    functor( Hdr, _, Arity ),
    Mtx = [Hdr|Rows],
    mtx_value_column_frequencies( Mtx, 1, Freqs ),
    write( freqs(Freqs) ), nl,
    findall( V-K, member(K-V,Freqs), Treqs ),
    sort( Treqs, Oreqs ),
    reverse( Oreqs, Rreqs ),
    findall( V-K, member(K-V,Rreqs), Preqs ),
    os_ext( pdf, Canc, PdfF ),
    os_path( pdfs, PdfF, PdfP ),
    % mlu_frequency_plot( Freqs, [interface(gg_bar),output(pdf(+PdfP))] ).
    atomic_list_concat( ['Number of patients per driver event for ',Clb,'.'], '', Main ),

    Mlu = [ interface(barplot), outputs([pdf(file= +PdfP, width=14, height=7)]),
            main= +Main, 
            xlab= "Drivers",
            ylab= "Patients",
            las=2, cex.names=0.5,
            post_call( <-( abline(h=Cut,col="#C23B22",lwd=2,lty=2)) )
          ],
    mlu_frequency_plot( Preqs, Mlu ),

    % aml_min00_muts_by_pnms.pdf :
    aml_add_mutations_column( Mtx, Stx ),
    aml_colour( airforceblue, AFBclr ),
    VLine = geom_segment(aes(x= -0.01, y = Cut, xend = Arity, yend = Cut), color= +AFBclr, linetype="dashed", size=0.4),
    % Vtext = annotate("text", x= -3, y= 60, hjust=0, vjust=0, label= "60"),
    % tmp <- 'data.frame'(tmpx=c(1,1),tmpy=c(2,3)),
    ( Canc == gbm -> 
        YTicks = scale_y_continuous(breaks = c(Cut,10,20,30,40,50,60))
        ;
        YTicks = scale_y_continuous(breaks = c(Cut,100,200,300,400))
    ),
    YTheme = theme(axis.text.y = element_text(color = c(+AFBclr, "black", "black", "black", "black"), size=c(11,9,9,9,9) ),
         axis.ticks.y = element_line(color = c(+AFBclr, "black", "black", "black", "black"),
                          size = c(1,1,1,1,1))),
    os_postfix( muts_by_pnms, PdfP, ByPnmsPdfP ),
    os_ext( pdf, ByPnmsPdfP, ByPnmsPdfS ),
    Pnms = [panel_theme(axes),gg_terms([VLine,YTicks,YTheme]),x_axis_colour(x_axis_aml)],
    gbn:gg_muts_by_pnms( Stx, ByPnmsPdfS, Pnms ),
        debug_call( etc, end, Canc ).
