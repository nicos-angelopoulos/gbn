
% :- lib( stoics_lib:io_lines/2 ).

/**  gbn_fix_dot( File )

Gobnilp doesnot quote headers of the form HLA-A
so the dots cannot be processed by dot executable.

INCOMPLETE. It was decided to just change the var names ... 
*/
gbn_fix_dot( File ) :-
	io_lines( File, Lines ),
	maplist( gbn_fix_dot_line, Lines, Fixed ),
	maplist( atom_codes, Atomed, Fixed ),
	tell( File ),
	maplist( writeln, Atomed ),
	told.

gbn_fix_dot_line( [], [] ).
gbn_fix_dot_line( [F|T], Fixed ) :-
	gbn_fix_dot_line( T, F, Fixed ).

gbn_fix_dot_line( [], Last, [Last] ).
gbn_fix_dot_line( [H|T], Prev, [FixPrev|More] ) :-
	gbn_fix_dot_code( H, Prev, T, FixPrev ),
    gbn_fix_dot_line( T, H, More ).
