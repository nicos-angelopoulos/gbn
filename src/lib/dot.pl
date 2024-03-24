
% :- lib(options).
% :- lib(stoics_lib:en_list/2).

dot_defaults(format(svg)).

/** dot( +DotsAndOpts )

Process a number dot files with dot -T png > $stem.png.

==
pupsh dot *rvr.tod
==

Now supports other outputs via: 
==
pupsh dot format=Fmt *rvr.tod
==

@author  nicos angelopoulos
@version 0.1 2020/07/15 specialised version for pack(gbn)

*/

dot( ArgS ) :-
     options_append( dot, ArgS, Opts, atoms(Files) ),
     options( format(Fmt), Opts ) ,
     debuc( dot, 'Format: ~w', [Fmt] ),
     maplist( dot_file(Fmt), Files ).

dot_file( Frm, File ) :-
     file_name_extension( Stem, _Dot, File ),
     file_name_extension( Stem, Frm, PngF ),
     atom_concat( '-T', Frm, TFrm ),
     debuc( dot, 'Doting file: ~w, onto: ~w', [File,PngF] ),
     AbsG = absolute_file_name( path(dot), _DotBin, [access(execute),file_errors(fail)] ),
     ( AbsG ->
         process_create( path(dot), [TFrm,'-o',PngF,File], [] )
         ;
         print_message( informational, gbn(exec_miss(dot,'cannot create net',File)) )
     ).
