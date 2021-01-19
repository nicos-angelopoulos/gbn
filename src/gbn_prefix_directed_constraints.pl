
% :- lib(stoics_lib:prefix_atom/2).

gbn_prefix_directed_constraints_defaults( bar_contra_self(true) ).

/** gbn_prefix_directed_constraints( +BnVars, +PfxPrs, +OutFile ).

For every pair of prefixes in PfxPrs:E-L, write all constrainsts that ban late BnVars (L prefixed)
from be the parents of early ones (E prefixed).

BnVars is a list of atomic Bn variables. OutFile is an open/3 compatible file for writing.

Opts 
  * bar_contra_self(Bar=true)
     whether to bar the early -> late edge for each BnVar. 
     These variables will invariably be connected in cases where each variable nearly always
     exists in either early or late stage

@author nicos angelopoulos
@version  0.1 2017/5/22

*/
gbn_prefix_directed_constraints( Vect, PfxPrsPrv, File, Args ) :-
    options_append( gbn_prefix_directed_constraints, Args, Opts ),
    options( bar_contra_self(Bar), Opts ),
    open( File, write, Out ),
    en_list( PfxPrsPrv, PfxPrs ),
    gbn_prefix_directed_constraints_stream( PfxPrs, Vect, Out, Bar ),
    close( Out ).

gbn_prefix_directed_constraints_stream( [], _Vect, _Out, _Bar ).
gbn_prefix_directed_constraints_stream( [Bef-Aft|T], Vect, Out, Bar ) :-
    include( prefix_atom(Bef), Vect, Befores ),
    include( prefix_atom(Aft), Vect, Afters ),
    maplist( gbn_directed_constraints_list_on_stream(Out,Befores), Afters ),
    gbn_prefix_directed_constraints_bar( Bar, Befores, Bef, Aft, Afters, Out ),
    gbn_prefix_directed_constraints_stream( T, Vect, Out, Bar ).

gbn_directed_constraints_list_on_stream( Out, Befores, After ) :-
    maplist( gbn_directed_constraint_on_stream(Out,After), Befores ).

gbn_directed_constraint_on_stream( Out, After, Before ) :-
    write( Out, '~' ),
    % maplist( replace_dot, [After,Before], [NoDotA,NoDotB] ),
    write( Out, Before ), write( Out, '<-' ),
    write( Out, After ), nl( Out ).

gbn_prefix_directed_constraints_bar( true, Befores, Bef, Aft, Afters, Out ) :-
    member( ABefore, Befores ),
    atom_concat( Bef, Stem, ABefore ),
    atom_concat( Aft, Stem, AnAfter ),
    memberchk( AnAfter, Afters ),
    gbn_directed_constraint_on_stream( Out, ABefore, AnAfter ),
    fail.
gbn_prefix_directed_constraints_bar( _, _Befores, _Bef, _Aft, _Afters, _Out ).
