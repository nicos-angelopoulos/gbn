
% :- lib(stoics_lib:n_digits_min/3).
% :- lib(stoics_lib:has_at_least/3).

gbn_mtx_mins_defaults( Defs ) :-
    Defs = [from(5),by(5),length(5)].

/** gbn_mtx_mins( +MtxF, -Files, +Opts ).

Create a number of matrix Files, curved out of input file MtxF.
A file is generate for each integer in Mins. The file in 
corresponding Min value is that which only contains columns
with at least Min number of 1s.

Files can be partially instantiated, to dictate the number of 
returned files.

Use debug(mtx_mins) to display the Min cutoffs, and the size of each matrix.

Opts
 * mins(Mins=[_,_,_,_,_])
   If ground, a list of the Mins to consider. 
   if non-ground, returns those that were generated 
   (can be a vars list to pass the implied Len)
 * from(From=5)
   starting point for Mins, when generating them
 * by(By=5)
   step
 * length(Len=5)
   if Mins is a variable, use this to determine length of Mins

==
?- 
==
@author nicos angelopoulos
@version  0:1 2019/3/24
@tbd debug() option in options_append/4 understands gbn(mtx_mins) ?

*/
gbn_mtx_mins( MtxIn, Files, Args ) :-
    options_append( gbn_mtx_mins, Args, Opts ),
    gbn_mtx_mins_gen( Files, Mins, Opts ),
    debug( gbn(mtx_mins), 'gbn mins: ~w', [Mins] ),
    mtx( MtxIn, Mtx, convert(true) ),
    mtx_lists( Mtx, Lists ),
    max_list( Mins, MaxMin ),
    number_codes( MaxMin, Codes ),
    length( Codes, Ndg ),
    maplist( gbn_mtx_lists_min(MtxIn,Ndg,Lists), Mins, Files ).

gbn_mtx_lists_min( InF, Ndg, Lists, Min, OutF ) :-
    n_digits_min( Ndg, Min, MinPad ),
    atom_concat( min, MinPad, Psfx ),
    os_postfix( Psfx, InF, OutF ),
    include( has_at_least(Min,1), Lists, MinLists ),
    mtx_lists( MinMtx, MinLists ),
    debug_call( gbn(mtx_mins), dims, Psfx/MinMtx ),
    mtx( OutF, MinMtx ).

gbn_mtx_mins_gen( Files, Mins, Opts ) :-
    options( mins(Mins), Opts ),
    ( ground(Mins) ->
        true
        ;
        ( var(Files) ->
            ( var(Mins) ->
                options( length(Len), Opts )
                ;
                length( Mins, Len )
            )
            ;
            length(Files,Len)
        ),
        options( from(From), Opts ),
        options( by(By), Opts ),
        findall( Min, ( between(1,Len,I),
                       Min is From + ((I-1) * By)
                     ), 
                        Mins )
    ).
