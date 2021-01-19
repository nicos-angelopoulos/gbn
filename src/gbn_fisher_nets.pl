
% :- lib( stoics_lib:en_list/2 ).
% :- lib( stoics_lib:map_list_options/4 ).

/** gbn_fisher_nets.
    gbn_fisher_nets(+Opts).

	Run a number of gbn_fisher_net/4 on files in current directory (version 2, as version 1 got lost).

	For now it assumes a single .dat file in current directory and runs this against all .bn files
	creating os_postfix/3 results files with postfix=fish.

	pupsh gbn_fisher_nets postfix=wb graph=bgcolor=white

	Options are passed to gbn_fisher_net/4 except for dir(Dir).

Opts 
 * dir(Dir='.')
   directory to change to

 * rec(Rec=false)
   allow recursive descend call

==
%  pupsh gbn_fisher_nets bground=white clr_theme=bic
==

Creates an undirected version of the above:
==
% pupsh gbn_fisher_nets bground=white clr_theme=bic type=graph postfix=biu
==

@author nicos angelopoulos
@version  0.1 2016/6/16

*/
gbn_fisher_nets :-
	gbn_fisher_nets( [] ).

gbn_fisher_nets( Args ) :-
	en_list( Args, Opts2 ),
	% maplist( string_fy_graph_colour, Opts1, Opts2 ),  % use backround/1 which auto strings it argument
	working_directory( Here, Here ),
	( select(dir(Dir),Opts2,Opts) -> working_directory(_Old,Dir) ; Opts = Opts2 ),
	os_files( Files ),
	include( os_ext(dat), Files, DatFs ),
	gbn_data_file( DatFs, DatF ),

	include( os_ext(bn), Files, BnFs ),
	map_list_options( gbn_fisher_net(DatF), BnFs, _Fishes, Opts ),  
            % do not use the gbn_fisher_net/3 version !
	working_directory( _, Here ).

gbn_data_file( [DatF], DatF ) :-
	!.
gbn_data_file( Other, _DatF ) :-
	throw( there_should_be_a_unique_data_file_instead_of:Other ).
