/** gbn_res_dir_dat_file( +Dir, -DatF ).

Locate a unique .dat file (DatF) within a given Dir-ectory.

Fails with a message if either none exists or more than one where located.

@author nicos angelopoulos
@version  0.1 2020/07/16
@tbd make sure all post-processors use this
@tbd control fail/error by flag ?

*/
gbn_res_dir_dat_file( Dir, DatF ) :-
    os_sel( os_files, ext(dat), Dats, dir(Dir) ),
	gbn_unique_dat_file( Dats, Dir, DatF ).

gbn_unique_dat_file( [Dat], _Dir, Dat ) :- !.
gbn_unique_dat_file( Dats, Dir, _ ) :-
	% throw( non_unique_dat_file(Dats) ).
    gbn_message( non_unique_dat(Dats,Dir) ),
    fail.
