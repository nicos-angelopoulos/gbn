
:- lib(mtx).
:- lib(real).
:- lib(b_real).
:- lib(options).

% :- debug(real).
:- lib(stoics_lib:kv_decompose/3).
:- <- library("RColorBrewer").

mtx_mut_hmap_test :-
    Mtx = [row(gene1,1,1,0,0,1,0,1,0,1,1),
           row(gene2,0,1,1,0,1,1,0,0,0,0)
          ],
    mtx_mut_hmap( Mtx, [] ).

mtx_mut_hmap_defaults( Defs ) :- 
    Types = [   as_mutational-boolean,
                hclust-oneof([clms,false]), rvar-atom,  rvar_rmv-boolean,
                % legend, atom or string % fixme: allow multitypes in type/n
                % height-number, width-number, dpi-integer,   % does options_append complain here ?
                % plot-var,
                stem-atom, x11-boolean
            ],
    Defs  = [ as_mutational(true),
              hclust(clms), 
              lbl_wt(background),
              lbl_mt(mutation),
              legend(bottom),
              legend_show(true),
              legend_font_size(10),
              % outputs: no default
              % plot(true), returns the plot term
              mtx_read(r),
              rvar(mtx_mut_hmap),
              rvar_rmv(true),
              stem(mtx_mut_hmap),
              x11(true),
              y_tick_size_x11(16),
              y_tick_size_cow(40),
              options_types(Types),
              col_mut("#CB181D"),
              col_wt("#08519C"),
              col_hmap([])
              ].

/** mtx_mut_hmap( +Mtx, +Opts ).

Display heatmap of a mutations matrix on screen or on file.

As of version 0.2 we support arbitrary integer value matrices. 
To do so, use AsMut == false and give a ClrH pairlist that will link each matrix value
to the color for drawing it.

This should have different backends, but currently only supports ggplot().

Opts 
  * as_mutational(AsMut=true)
    whether the matrix contains mutational data (when false arbitary integer values are assumed)
  * col_hmap(ClrH=[])
    list of value-"Colour" for each value in mtx (can contain non mtx values).
    If this is given it overrides ClrM and ClrW.
  * col_mut(ClrM="#CB181D")),
    colour from mutations (mtx values = 1)
  * col_wt(ClrW="#08519C")),
    colour from background (mtx values = 0)
  * dpi(Dpi=300)
    dpi for output files
  * hclust(Hclust=clms)
    whether and what to cluster, currently either clms or false
  * height(Height)
    number for Height in mm for file plotting, defaults to: 10 + (10 * Nr) + LYpad,
  * legend(LegPos=bottom)
    atom or string of legend position (set to false/FALSE for no legend)
  * legend_show(ShowLeg=true)
    whether to show the legend
  * legend_font_size(LegFntSz=10)
    font size for the legend text
  * lbl_mt(Lmt=mutation)
    label for 1 values (also expands to LmtX for values X > 1)
  * lbl_wt(Lfg=background)
    label for 0 values
  * mtx(Mtx)
    the input Matrix for the mtx_mut_hmap/1 version.<br>
    Has no effect on mtx_mut_hmap/2 version. Alternative value: _pl_.
  * mtx_read(Mfc=r)
    interface for reading-in the data csv, if Mtx is a file.<br>
    The default is _true_ but is ignored if Mtx is not a file.
  * outputs(OutS)
    an output or list of outputs, each being either atomic or compound. 
    The functor, or atom, of each output should be an file extension understood by R's ggsave().
    Arguments of the compound should be = pairs also understood by ggsave().
    Default values for height, width and dpi (300) are given. In the case of height and width
    these are calculated as a function of the input matrix dimensions. Eg. OutS = png(height=200)
    produces a png image with non-default height at 200 mm.
  * plot(Plot)
    returns the plot term. if present, RvRmv is ignored and data frame is kept (as it will be needed when plotting Plot)
  * rvar(Rvar=mtx_mut_hmap)
    variable to use for intermediate result
  * rvar_rmv(RvRmv=true)
    removes R var at end of call, see options_rvar_rmv/2 in pack b_real
  * stem(Stem=mtx_mut_hmap)
    set to stem for output files  
  * width(Width)
    number for Width in mm for file plotting, defaults to: 20 + Nc + LXpad,
  * x11(X11=true)
    whether to plot on screen
  * y_tick_size_cow(YtC=40)
     size for y_ticks on the cowplots
  * y_tick_size_x11(YtX=16)
     size for y_ticks on X11 (and individual heatmaps)

@author  nicos angelopoulos
@version 0.1 2017/02/14
@version 0.2 2022/03/19, support non-mutational hmaps
@version 0.3 2024/05/08, reduce PL <- R interactions around hclust()
@version 0.4 2024/08/17, single arg version
@tbd     make default rvar a non-existant variable

*/
mtx_mut_hmap( Args ) :-
    options_append( mtx_mut_hmap, Args, Opts, pack(sanger) ),
    options( mtx(Mtx), Opts ),
    mtx_mut_hmap_opts( Mtx, Opts ).
     
mtx_mut_hmap( Mtx, Args ) :-
    options_append( mtx_mut_hmap, Args, Opts, pack(sanger) ),
    mtx_mut_hmap_opts( Mtx, Opts ).

mtx_mut_hmap_opts( MtxIn, Opts ) :-
    options( as_mutational(Bin), Opts ),
    options( rvar(Rvar), Opts ),
    options( hclust(Hcl), Opts ),
    options( legend_show(LegShow), Opts ),
    upcase_atom( LegShow, ShowLeg ),
    options( legend_font_size(LegFntSz), Opts ),
    options( mtx_read(Mfc), Opts ),
    mtx_mut_hmap_mtx( Mfc, MtxIn, Rvar ),

    mtx_mut_hmap_cluster( Hcl, Rvar ),
    % mtx_mut_hmap_cluster( false, Rvar ),
    Nr <- nrow(Rvar),
    Nc <- ncol(Rvar),
    Rwms <- rownames(Rvar),
    atomic_list_concat( [Rvar,df], '_', Df ),
    Df <- data.frame( x=integer(), y=character(), m=character(), stringsAsFactors='FALSE' ),
    % fixme: convert this to a recursion
    % there is a version in $local@ampelos, but is not correct
    options( lbl_wt(Lwt), Opts ),
    options( lbl_mt(Lmt), Opts ),
    mtx_mut_hmap_df_rows( Rwms, 1, Nc, Bin/Lwt/Lmt, Rvar, Df ),

    % trace  % fixme: when x11 is off the plot should be passed to ggsave()
    % <- print( Df ),
    Rvar <- Df,
    Rvar$y <- as.character(Rvar$y),
    Rvar$y <- factor(Rvar$y, levels=rev(Rwms)),
    r_remove( Df ),
    options( legend(LegPos), Opts ),
    mtx_mut_hmap_leg_pads( LegPos, LXpad, LYpad, ActLegPos ),
    PreLvls <- levels(as.factor(Rvar$m)),
    ( catch(maplist(atom_number,PreLvls,LvlsN),_,fail) ->
          % if all values are numeric, then sort them by numeric value and not lexigocraphically
          sort( LvlsN, LvlsO ),
          maplist( atom_number, NewLvls, LvlsO ),
          Rvar$m <- factor( as.factor(Rvar$m), levels = NewLvls )
          ;
          % else leave them as they are
          true
    ),
    LvlsPrv <- levels(as.factor(Rvar$m)),
    (is_list(LvlsPrv) -> LvlsPrv = Lvls; Lvls = [LvlsPrv]), % fixme: throw error if not a list...
    options( col_hmap(ClrH), Opts ),
    mtx_mut_hmap_colours( ClrH, Lvls, ClrsL, Opts ),
    Clrs =.. [c|ClrsL],   % because cowplot needs them un-magicked...
    Gp = ggplot(Rvar) + geom_tile( aes(x=x,y=y,fill=m), 'show.legend'=ShowLeg ) 
           + scale_fill_manual( values=Clrs)
           + theme( legend.position= +ActLegPos
                    , axis.text = element_text(size = TickSize)
                    , axis.title.x=element_blank()
                    , axis.text.x=element_blank()
                    , axis.ticks.x=element_blank()
                    , axis.title.y=element_blank()
                    , legend.title=element_blank()    % removes legend title
                    , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                    , panel.background = element_blank()
                    , panel.spacing= unit(2,"cm")
                  )
          + guides(fill = guide_legend(order = 2, keywidth = 0.8, keyheight = 1.6, 
                                label.theme = element_text(size = LegFntSz,face = "italic",angle = 0)
                                      )
                  ),
    options( [y_tick_size_x11(YtX),y_tick_size_cow(YtC)], Opts ),
    findall( Gp, member(TickSize,[YtX,YtC]), [GpX,GpC] ),
    options_call( x11(true), real:(<-(print(GpX))), Opts ),
    Height is (10 * (Nr + 1)) + LYpad,
    Width  is min( max( 20 + (Nc/4) + LXpad, 70 ), 1024 ),
    append( Opts, [height(Height),width(Width),dpi(300)], SaveOpts ),
    mtx_mut_hmap_save( GpX, SaveOpts ),
    % nicos.fixme.delete.me.24.08.16
    ( memberchk(plot(Plot),Opts) ->
        Plot = GpC
        ;
        options_rvar_rmv( Rvar, Opts )
    ).

mtx_mut_hmap_mtx( r, MtxIn, Rvar ) :-
    atomic( MtxIn ),
    !,
    Rvar <- read.csv( +MtxIn, header='FALSE', row.names=1 ).
mtx_mut_hmap_mtx( _, MtxIn, Rvar ) :-
    mtx( MtxIn, Mtx ),
    findall( Rwn-DtRow, ( member(Row,Mtx), 
                          Row =.. [Rfun,Rwn|Dt],
                          DtRow =.. [Rfun|Dt]
                        ),
                            RDs ),
    kv_decompose( RDs, Rwns, Rtx ),
    Rvar <- Rtx,
    rownames(Rvar) <- Rwns.

mtx_mut_hmap_df_rows( [], _I, _Nc, _Bin, _Rv, _Df ).
mtx_mut_hmap_df_rows( [Rwm|Rwms], I, Nc, Bin, Rvar, Df ) :-
    /* was
    findall( _,         (   nth1(Rn,Rwms,Rwm),
                            between(1,Nc,Cn),
                            Ri is ((Rn - 1) * Nc ) + Cn,
                            Df[Ri,1] <- Cn,
                            Df[Ri,2] <- +Rwm,
                            Mut <- Rvar[Rn,Cn],
                            mtx_mut_hmap_clr_value( Bin, Mut, Mfc ),
                            Df[Ri,3] <- Mfc
                         ),
                            _ ),
                            */
     mtx_mut_hmap_df_vals( 1, Nc, Rwm, I, Bin, Rvar, Df ),
     J is I + 1,
     mtx_mut_hmap_df_rows( Rwms, J, Nc, Bin, Rvar, Df ).

mtx_mut_hmap_df_vals( I, Nc, _Rwm, _Rn, _Bin, _Lwt, _Lmt,  _Rv, _Df ) :-
     Nc < I,
     !.
mtx_mut_hmap_df_vals( Cn, Nc, Rwm, Rn, Bin, Lwt, Lmt, Rvar, Df ) :-
     Ri is ((Rn - 1) * Nc ) + Cn,
     Df[Ri,1] <- Cn,
     Df[Ri,2] <- +Rwm,
     Mut <- Rvar[Rn,Cn],
     mtx_mut_hmap_clr_value( Bin, Lwt, Lmt, Mut, Mfc ),
     Df[Ri,3] <- Mfc,
     Co is Cn + 1,
     mtx_mut_hmap_df_vals( Co, Nc, Rwm, Rn, Bin, Lwt, Lmt, Rvar, Df ).

mtx_mut_hmap_colours( [H|T], Lvls, Clrs, _Opts ) :-
     !,
     Pairs = [H|T],
     findall(+Clr, (member(Lvl,Lvls),(memberchk(Lvl-Clr,Pairs)->true;throw(missing_colour_spec(Lvl)))), Clrs ),
     write( colours(Clrs) ), nl.
mtx_mut_hmap_colours( _, Lvls, Clrs, Opts ) :-
    length( Lvls, LvlsLen ),
    RedAtms <- brewer.pal(9,"Reds"),
    maplist( atom_string, RedAtms, Reds ),
    options( [col_mut(ClrM),col_wt(ClrW)], Opts ),
    ( memberchk(background,Lvls) ->
        ( LvlsLen > 10 -> throw( too_many_levels_in_discrete_variable_for_mut_hmap(Lvls) )
                          ;
                          ( LvlsLen =:= 2 ->
                                Clrs = [ClrW,ClrM]
                                ;
                                Lim is LvlsLen - 1,
                                length( RedClrs, Lim ),
                                once( append(_,RedClrs,Reds) ),
                                Clrs = [ClrW|RedClrs]
                          )
        )
        ;
        ( LvlsLen > 9 -> throw( too_many_levels_in_discrete_variable_for_mut_hmap(Lvls) )
                         ; 
                         length( Clrs, LvlsLen ),
                         once( append(_,Clrs,Reds) )
        )
    ).

mtx_mut_hmap_clr_value( false, _Lwt, _Lmt, Mut, Mfc ) :-
     Mfc = Mut.
mtx_mut_hmap_clr_value( true, Lwt, Lmt, Mut, +(Mfc) ) :-
     %  mutation = 1, background 0
     %  however, we allow for other values as there might be other discrete values in the dataset
     %  thus we interpret 0 as background and mutI for each value > 0
     ( Mut =:= 0 ->   % 22.03.19 this used to be Mut < 1, which is probably more correct if we
                      % we insist that this a mutational matrix (AsMut==true), even if there are negative
                      % numbers in the matrix. 
                      % This impelentation is "better" in highlighting something is gone wrong, so it is left in
          Mfc = Lwt
          ;
          ( Mut =:= 1 -> 
               Mfc = Lmt
               ;
               atom_concat( Lmt, Mut, Mfc )
          )
    ).

mtx_mut_hmap_save( Ggl, Opts ) :-
    memberchk( outputs(OutS), Opts ),
    !,
    options( stem(Stem), Opts ),
    options( [height(H),width(W),dpi(D)], Opts ),
    en_list( OutS, Outs ),
    maplist( mtx_mut_hmap_save_on(Ggl,Stem,[height=H,width=W,dpi=D]), Outs ).
mtx_mut_hmap_save( _Opts ).

mtx_mut_hmap_save_on( Ggl, Stem, Defs, Out ) :-
    atomic( Out ),
    !, 
    mtx_mut_hmap_save_on_opts( Ggl, Stem, Out, Defs ).
mtx_mut_hmap_save_on( Ggl, Stem, Defs, Out ) :-
    Out =.. [Func|Args],
    append( Args, Defs, Opts ),
    mtx_mut_hmap_save_on_opts( Ggl, Stem, Func, Opts ).

mtx_mut_hmap_save_on_opts( Ggl, Stem, Ext, Opts ) :-
    os_ext( Ext, Stem, File ),
    memberchk( height=Height, Opts ),
    memberchk( width=Width, Opts ),
    memberchk( dpi=Dpi, Opts ),
    <- ggsave( +File, plot=Ggl, height=Height, width=Width, units="mm", dpi=Dpi, limitsize='FALSE' ).

mtx_mut_hmap_leg_pads( LegPos, LXpad, LYpad, ActLegPos ) :-
    ( string(LegPos) -> atom_string( AtmLegPos, LegPos ) ; AtmLegPos = LegPos ),
    downcase_atom( AtmLegPos, LowLegPos ),
    mtx_mut_hmap_leg_pads_atom( LowLegPos, LXpad, LYpad, ActLegPos ).

mtx_mut_hmap_leg_pads_atom(false,  0,  0,  false).
mtx_mut_hmap_leg_pads_atom(bottom, 0,  8, bottom).
mtx_mut_hmap_leg_pads_atom(top,    0,  8, bottom).
mtx_mut_hmap_leg_pads_atom(left,  30,  0,   left).
mtx_mut_hmap_leg_pads_atom(right, 30,  0,  right).
mtx_mut_hmap_leg_pads_atom(true, X, Y, Act) :-
    mtx_mut_hmap_leg_pads_atom( bottom, X, Y, Act ).

mtx_mut_hmap_cluster( clms, Rvar ) :-
    rv_tmp_dist <- dist(t(Rvar)),
    rv_tmp_ord <- hclust(rv_tmp_dist)$order,
    % fixme: the following line i think was responsible for "Killed" crash, when called from a findall/3 call:
    % Ord <- hclust( dist(t(Rvar)) )$order,
    Rvar <- Rvar[*,rv_tmp_ord],
    <- remove(rv_tmp_dist),
    <- remove(rv_tmp_ord).
mtx_mut_hmap_cluster( false, _Rvar ).
