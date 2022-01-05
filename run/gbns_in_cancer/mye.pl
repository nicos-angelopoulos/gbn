
:- use_module(library(lib)).

:- lib(gbn).
:- lib(debug_call).    % debuc/1, debu/3.

:- debuc(mye).

/** mye.

Run myeloma experiments as reported on paper.

Data copied from myeloma/rea.18/wgs/data/rea_min20.dat.

*/

mye :-
        debuc( mye, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/mye/mye_min20.dat'),
    E = 2,
    GBNOpts = [ data(DatF),
                setting(edge_penalty,E),
                dir(OsOdir)
                % debug(true)
              ],
    gbn( GBNOpts ),
        debuc( mye, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debuc( mye, start, fam_hmaps ),
    gbn_fam_hmaps( [dir(OsOdir),not(true)] ),   % only for binaries
        debuc( mye, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debuc( mye, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debuc( mye, end, true ).   

display_var_as( hyper, "HRD" ).
display_var_as( t4_14, "t(4;14)" ).
display_var_as( t11_14, "t(11;14)" ).
display_var_as( t14_16, "t(14;16)" ).

/** mye_full.

Run GBN on the myeloma cancer dataset for \epsilon =1...10.

*/
mye_full :-
    
        debuc( mye, start, true ),
    DatF = pack('gbn/data/gbns_in_cancer/mye/mye_min20.dat'),
    numlist( 1, 10, Es ),
    GBNOpts = [ data(DatF),
                multiple(edge_penalty,Es),
                dir(OsOdir)
              ],
    gbn( GBNOpts ),
        debuc( mye, start, fisher_nets ),
    gbn_fisher_nets( dir(OsOdir) ),    
        debuc( mye, start, fam_hmaps ),
    gbn_fam_hmaps( [dir(OsOdir),not(true)] ),   % only for binaries
        debuc( mye, start, gates_nets ),
    gbn_gates_nets( dir(OsOdir) ),
        debuc( mye, start, svg_legend ),
    gbn_svg_legend( dir(OsOdir) ),
        debuc( mye, end, true ).   
