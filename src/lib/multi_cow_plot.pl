
:- lib(r(cowplot)).   
        % fixme: make this to a soft load ? if so make sure you deal with missing in the code
% :- lib(stoics_lib:letter_strings/3).

multi_cow_plot_defaults( Defs ) :-
        Defs = [   ext(png),
            % hjust(-0.01), % default is -0.5
            aspect(0.8),
            hjust(0.10),   % larger values move labels to the left 
            labels(upper),
            ncol(2),
            stem(multi_cow_plot),
            vjust(0.9),
            % vjust(0.5)   % smaller means higher up, use this with scale 0.96
            height(_)
            ].

/** multi_cow_plot( +Plvs, +Opts ).

Create a single multi plot of many ggplot R variables (Plvs) by using library(cowplot).

Opts 
  * aspect(Aspect=1.8)
     value of the base_aspect_ratio parameter of cow plot

  * ext(Ext=png)
     type of file to save on

  * height(Height=H)
     height of the plot, if variable returns the height used
     (H is max( ( ( (Lenvs + 1) // 2 ) * 1.5) * 0.0393701, 30) )
     
  * hjust(JustH= -1)
     horizontal justification (plot_grid() param.)

  * labels(Lbls=upper)
     also can use lower, true (=upper), auto  and false. Else pass a list wihcih is passed to the plotter
     auto passes "auto" as labels which is recognised by the draing function (plot_grid())

  * ncol(Ncol=2)
     number of columns (plot_grid() param.)

  * stem(Stem=multi_cow_plot)
     stem of file to save on

  * vjust(JustV= -1)
     vertical justification (plot_grid() param.)

@author nicos angelopoulos
@version  0.1 2017/02/15
@version  0.2 2017/07/18,  height/1 and apect/1 options, untested...
@tbd  maximum per page
@tbd  exclude labels
@tbd  rvar-ise mplot

*/
multi_cow_plot( Plvs, Args ) :-
    options_append( multi_cow_plot, Args, Opts ),
    length( Plvs, Lenvs ),
    options( labels(LblsTkn), Opts ),
    multi_cow_plot_labels( LblsTkn, Lenvs, Lbls ),
    options( [ncol(Ncol),hjust(JustH),vjust(JustV)], Opts ),
    % Rargs = [scale=0.96,labels=Lbls,ncol=Ncol,hjust=JustH,vjust=JustV], % use scale= to fit labels in margins
    Rargs = [labels=Lbls,ncol=Ncol,hjust=JustH,vjust=JustV,align="v"],
    append( Plvs, Rargs, PgArgs ),
    Grid =.. [plot_grid|PgArgs],
    % <- theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")),  maybe theme_cow ? 
    mplot <- Grid,
    'mplot$layout$clip[mplot$layout$name == "panel"]' <- "off",
    options( [stem(Stem),ext(Ext)], Opts ),
    os_ext( Ext, Stem, File ),
    options( aspect(Aspect), Opts ),
    options( height(H), Opts ),
    ( var(H) -> H is max( ( ( (Lenvs + 1) // 2 ) * 1.5) * 0.0393701, 30)
              ; true ),
    % W is max( 20 + (Nc/4) + LXpad, 70 ),
    % <- save_plot( +File, mplot, ncol=2, base_height=H, base_width=W ),
    <- save_plot( +File, mplot, ncol=2, base_height=H, base_aspect_ratio=Aspect, limitsize='FALSE' ),
    r_remove( mplot ).

multi_cow_plot_labels( LblsTkn, Lnvs, Lbls ) :-
    multi_com_plot_labels_expand_known(LblsTkn,Lnvs,Lbls),
    !.
multi_cow_plot_labels( LblsTkn, Lnvs, Lbls ) :-
    Use = lower, % defaulty
    gbn_message( cow_plot_lbl(LblsTkn,Use) ), 
    multi_com_plot_labels_expand_known(Use,Lnvs,Lbls).

multi_com_plot_labels_expand_known( auto, _, Auto ) :- atom_string( auto, Auto ).  % fixme: checkme
multi_com_plot_labels_expand_known( false, _,  [] ). % fixme: checkme
multi_com_plot_labels_expand_known( true, Lenvs, Lbls ) :-
    multi_cow_plot_labels( upper, Lenvs, Lbls ).
multi_com_plot_labels_expand_known( upper, Lenvs, Letts ) :-
    letter_strings( "A", Lenvs, Letts ).
multi_com_plot_labels_expand_known( lower, Lenvs, Letts ) :-
    letter_strings( "a", Lenvs, Letts ).
