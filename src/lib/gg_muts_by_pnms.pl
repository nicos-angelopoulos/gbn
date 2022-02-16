
:- lib(mtx).
:- lib(real).
:- lib(b_real).     % vectors_subed_gg_bar_plot/4, c25/1.
:- lib(os_lib).     % os_dir_stem_ext/3
:- lib(options).

gg_muts_by_pnms_defaults( [ odir('.'), ext(pdf),
                            % x_axis_colour(hmrn:co_blood),
                            x_axis_colour(false),
                            extra_legend_position(c(88,102,20,50)),
                            plot_save_width(21),
                            plot_title(Main),
                            plot_leg_title(LegTitle)
                     ] ) :-
    muts_by_pnms_titles( Main, LegTitle ).

/** gg_muts_by_pnms( +Mtx, +Stem, -Opts ).

   For matrix Mtx ggplot histogram of each mutational 
   variable coloured by mutations column. 

    * ext(Ext=pdf)
      extension for output file (see os_dir_stem_ext/2).

    * extra_legend_position(Lp=c(88,102,20,50))
      position for the extra legend (if one is required)

    * odir(Odir='')
      output directory (as per os_mill/4)

    * plot_title(Main='Number of mutations for each driver event.'),
       title for plot (appears at the top)

    * plot_leg_title(LegTitle='Patient\nmutations')
       title for legend

    * plot_save_width(Width=21)
       width in centimeters for the saved output

    * x_axis_colour(XaClrs=false)
      colours for the ticks on x-axis


Options are passed to header_cnm_muts/3 and vectors_subed_gg_bar_plot/4 (and via that to gg_bar_plot/2).

Rows with 0 mutations are now removed (20.7.29) previously a zero would appear on the
legend if 0 mutation rows existed in the dataset.

@author nicos angelopoulos
@version  0:2 2020/08.29
@see header_cnm_muts/3
@see vectors_subed_gg_bar_plot/4

*/
gg_muts_by_pnms( MtxIn, Stem, Args ) :-
    options_append( gg_muts_by_pnms, Args, Opts ),
    mtx( MtxIn, MtxAll ),
    MtxAll = [Hdr|Rows],
    Hdr =.. [_|Cnms],
    % include( cnm_symbol, Cnms, Symbs ),
    header_cnm_muts( Hdr, ClassCnm, Opts ),

    % fixme: this is a temporsCCary fix, need to have option->hook that defines Symbs
    %        Symbs are the columns that hold mutations/events that will be plotted
    once( nth1(PosCC,Cnms,ClassCnm,Symbs) ), 
    exclude( gg_muts_zero_summed(PosCC), Rows, NZows ),
    Mtx = [Hdr|NZows],
    c25( C25 ),
    options( extra_legend_position(c(P1,P2,P3,P4)), Opts ),
    % extra_legend_position(88,102,20,50)
    options( x_axis_colour(XaClr), Opts ),
    options( plot_title(Main), Opts ),
    options( plot_leg_title(LegTitle), Opts ),
    GGopts = [sort_x(totals),mtx(Mtx),x_axis_colour(XaClr),
              % fixme Patients -> samples
              labels("\nDrivers","Patients",Main),
              colours(C25),legend_title(LegTitle),
              extra_legend_position(P1,P2,P3,P4)
            | Opts
            ],
    user:vectors_subed_gg_bar_plot( Symbs, 1, ClassCnm, GGopts ),
    os_dir_stem_ext( File, [stem(Stem)|Opts] ),
     options( plot_save_width(Width), Opts ),
    <- ggsave( +File, width=Width ).

gg_muts_zero_summed( Pos, Row ) :-
    arg( Pos, Row, 0 ).

muts_by_pnms_titles( Main, Leg ) :-
    Leg  = "Patient\nmutations",
    % Main = "Total number of mutations for drivers coloured by patient mutations total\n".
    Main ='Number of mutations for each driver event.'.

% cnm_symbol( A, A ).

/** header_cnm_muts( +Hdr, -Cnm, +Opts ).

Identify the mutations column in Opts or Header. If none is found 
then 'Mutations', 'Patient_mutations' and 'Patient mutations' are 
looked for in Hdr.

First 

Opts 
 * cnm_muts(Cnm)
   user provided column name for mutations column (sum of mutations for each patient)

*/
header_cnm_muts( _Hdr, Cnm, Args ) :-
    en_list( Args, Opts ),
    memberchk( cnm_muts(Cnm), Opts ),
    !.
header_cnm_muts( Hdr, Cnm, _Opts ) :-
    cnm_muts_candidate( Cnm ),
    arg( _, Hdr, Cnm ),
    !.
header_cnm_muts( Hdr, _Cnm, Opts ) :-
    throw( cannot_locate_mutations_column_in_header_with_options(Hdr,Opts) ).

cnm_muts_candidate( 'Mutations' ).
cnm_muts_candidate( 'Patient_mutations' ).
cnm_muts_candidate( 'Patient mutations' ).
