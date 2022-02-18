
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
    Types = [   hclust-oneof([clms,false]), rvar-atom,  rvar_rmv-boolean,
                % legend, atom or string % fixme: allow multitypes in type/n
                % height-number, width-number, dpi-integer,   % does options_append complain here ?
                % plot-var,
                stem-atom, x11-boolean
            ],
    Defs  = [ hclust(clms), 
              legend(bottom),
              legend_show(true),
              legend_font_size(10),
              % outputs: no default
              % plot(true), returns the plot term
              rvar(mtx_mut_hmap),
              rvar_rmv(true),
              stem(mtx_mut_hmap),
              x11(true),
              y_tick_size_x11(16),
              y_tick_size_cow(40),
              options_types(Types)
              ].

/** mtx_mut_hmap( +Mtx, +Opts ).

Display heatmap of a mutations matrix on screen or on file.

This should have different backends, but currently only supports ggplot().

Opts 
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
     size for y_ticks on the cow plots
  * y_tick_size_x11(YtX=16)
     size for y_ticks on X11 (and individual heatmaps)

@author nicos angelopoulos
@version  0.1 2017/02/14
@tbd    make default rvar a non-existant variable

*/
mtx_mut_hmap( MtxIn, Args ) :-
    options_append( mtx_mut_hmap, Args, Opts, pack(sanger) ),
    options( rvar(Rvar), Opts ),
    options( hclust(Hcl), Opts ),
    options( legend_show(LegShow), Opts ),
    upcase_atom( LegShow, ShowLeg ),
    options( legend_font_size(LegFntSz), Opts ),
    mtx( MtxIn, Mtx ),
    % mtx( MtxIn, MtxRev ),
    % reverse( MtxRev, Mtx ),
    findall( Rwn-DtRow, ( member(Row,Mtx), 
                          Row =.. [Rfun,Rwn|Dt],
                          DtRow =.. [Rfun|Dt]
                        ),
                            RDs ),
    kv_decompose( RDs, Rwns, Rtx ),
    Rvar <- Rtx,
    rownames(Rvar) <- Rwns,
    mtx_mut_hmap_cluster( Hcl, Rvar ),
    Nr <- nrow(Rvar),
    Nc <- ncol(Rvar),
    Rwms <- rownames(Rvar),
    atomic_list_concat( [Rvar,df], '_', Df ),
    Df <- data.frame( x=integer(), y=character(), m=character(), stringsAsFactors='FALSE' ),
    findall( _,         ( % between(1,Nr,Rn),
                            nth1(Rn,Rwms,Rwm),
                            between(1,Nc,Cn),
                            Ri is ((Rn - 1) * Nc ) + Cn,
                            Df[Ri,1] <- Cn,
                            % atom_string( Rwm, Swm ),
                            Df[Ri,2] <- +Rwm,
                            Mut <- Rvar[Rn,Cn],
                            % 17.10.04: select reds for mutations for that 1
                            %           make sure though in colours that there is no blue if thee is no zero

                            %  mutation = 1, background 0
                            %  however, we allow for other values as there might be other discrete values in the dataset
                            %  thus we interpret 0 as background and mutI for each value > 0
                            ( Mut > 0 -> 
                                        ( Mut =:= 1 -> 
                                            Mfc = mutation
                                            ;
                                            atom_concat( mutation, Mut, Mfc )
                                        )
                                       ; 
                                       Mfc = background 
                            ),
                            % ( Mut =:= 1 -> Mfc = mutation; Mfc = background ),  % if labels are changed make sure colours coordinate
                            Df[Ri,3] <- + Mfc
                         ),
                            _ ),
    % Df$y <- as.factor(Df$y),»
    % <- print( unique(Df$y) ),
    % Df$y <- as.character(Df$y),
    % Df$y <- factor(Df$y, levels=unique(Df$y)),
    % Df$m <- as.factor(Df$m),
    Rvar <- Df,
    Rvar$y <- as.character(Rvar$y),
    % Rvar$y <- factor(Rvar$y, levels=unique(Rwns)),
    Rvar$y <- factor(Rvar$y, levels=rev(Rwns)),
    % Rvar$y <- factor(Rvar$y, levels=unique(Rvar$y)),
    % <- Rvar$y,
    r_remove( Df ),
    options( legend(LegPos), Opts ),
    mtx_mut_hmap_leg_pads( LegPos, LXpad, LYpad, ActLegPos ),
    GotMuns <- Rvar$m,
    sort( GotMuns, GotMprv ),
    ( is_list(GotMprv) -> GotMprv = GotM; GotM = [GotMprv] ),
    length( GotM, GotMlen ),
    RedAtms <- brewer.pal(9,"Reds"),
    maplist( atom_string, RedAtms, Reds ),
    ( memberchk(background,GotM) ->
        ( GotMlen > 10 -> throw( too_many_levels_in_discrete_variable_for_mut_hmap(GotM) )
                          ;
                          ( GotMlen =:= 2 ->
                                ClrsL = ["#08519C","#CB181D"]
                                ;
                                Lim is GotMlen - 1,
                                length( RedClrs, Lim ),
                                once( append(_,RedClrs,Reds) ),
                                ClrsL = ["#08519C"|RedClrs]
                          )
        )
        ;
        ( GotMlen > 9 -> throw( too_many_levels_in_discrete_variable_for_mut_hmap(GotM) )
                         ; 
                         length( ClrsL, GotMlen ),
                         once( append(_,ClrsL,Reds) )
        )
    ),
    Clrs =.. [c|ClrsL],   % because cows plot needs them un-magicked...
    Gp = ggplot(Rvar) + geom_tile( aes(x=x,y=y,fill=m), 'show.legend'=ShowLeg ) 
           % + scale_fill_manual(values=c("#CB181D","#08519C"))
           + scale_fill_manual( values=Clrs )   % values=c("#08519C","#CB181D","#000000") % CAREFULL order depends on lex order of $m
           + theme( legend.position= +ActLegPos
                    , axis.text = element_text(size = TickSize)
                    , axis.title.x=element_blank()
                    , axis.text.x=element_blank()
                    , axis.ticks.x=element_blank()
                    , axis.title.y=element_blank()
                    , legend.title=element_blank()    % removes legend title
                    , panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                    , panel.background = element_blank()
                    , panel.margin= unit(c(2,2,2,2), "cm")
                  )
          + guides(fill = guide_legend(order = 2, keywidth = 0.8, keyheight = 1.6, 
                                       label.theme = element_text(size = LegFntSz,face = "italic",angle = 0)  % colour="red", 
                                      )
                   % fill = guide_legend(order = 2, override.aes = list(shape = 21, size=3)) , shape = guide_legend( order= 1)
                  )
        ,
    options( [y_tick_size_x11(YtX),y_tick_size_cow(YtC)], Opts ),
    findall( Gp, member(TickSize,[YtX,YtC]), [GpX,GpC] ),

    options_call( x11(true), real:(<-(print(GpX))), Opts ),

    % Height is max( (10 * (Nr + 1)) + LYpad, 1200 ),
    Height is (10 * (Nr + 1)) + LYpad,
    Width  is max( 20 + (Nc/4) + LXpad, 70 ),
    append( Opts, [height(Height),width(Width),dpi(300)], SaveOpts ),
    mtx_mut_hmap_save( SaveOpts ),
    ( memberchk(plot(Plot),Opts) ->
        Plot = GpC
        ;
        options_rvar_rmv( Rvar, Opts )
    ).

mtx_mut_hmap_save( Opts ) :-
    memberchk( outputs(OutS), Opts ),
    !,
    options( stem(Stem), Opts ),
    options( [height(H),width(W),dpi(D)], Opts ),
    en_list( OutS, Outs ),
    maplist( mtx_mut_hmap_save_on(Stem,[height=H,width=W,dpi=D]), Outs ).
mtx_mut_hmap_save( _Opts ).

mtx_mut_hmap_save_on( Stem, Defs, Out ) :-
    atomic( Out ),
    !, 
    mtx_mut_hmap_save_on_opts( Stem, Out, Defs ).
mtx_mut_hmap_save_on( Stem, Defs, Out ) :-
    Out =.. [Func|Args],
    append( Args, Defs, Opts ),
    mtx_mut_hmap_save_on_opts( Stem, Func, Opts ).

mtx_mut_hmap_save_on_opts( Stem, Ext, Opts ) :-
    os_ext( Ext, Stem, File ),
    memberchk( height=Height, Opts ),
    memberchk( width=Width, Opts ),
    memberchk( dpi=Dpi, Opts ),
    <- ggsave( +File, height=Height, width=Width, units="mm", dpi=Dpi ).

mtx_mut_hmap_leg_pads( LegPos, LXpad, LYpad, ActLegPos ) :-
    ( string(LegPos) -> atom_string( AtmLegPos, LegPos ) ; AtmLegPos = LegPos ),
    downcase_atom( AtmLegPos, LowLegPos ),
    mtx_mut_hmap_leg_pads_atom( LowLegPos, LXpad, LYpad, ActLegPos ).

mtx_mut_hmap_leg_pads_atom( false,  0,  0,  false ).
mtx_mut_hmap_leg_pads_atom( bottom, 0,  8, bottom ).
mtx_mut_hmap_leg_pads_atom( top,    0,  8, bottom ).
mtx_mut_hmap_leg_pads_atom( left,  30,  0,   left ).
mtx_mut_hmap_leg_pads_atom( right, 30,  0,  right ).
mtx_mut_hmap_leg_pads_atom( true, X, Y, Act ) :-
    mtx_mut_hmap_leg_pads_atom( bottom, X, Y, Act ).

mtx_mut_hmap_cluster( clms, Rvar ) :-
    Ord <- hclust( dist(t(Rvar)) )$order,
    Rvar <- Rvar[*,Ord].
mtx_mut_hmap_cluster( false, _Rvar ).
