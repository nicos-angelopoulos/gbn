
:- lib(dot_read/2).
:- lib(dot_write/2).

/** dot_fix( +DotF ).

Simply reading the DotF and re-writing fixes 
some of the incompatibilities.

*/
dot_fix( DotF ) :-
	dot_read( DotF, Graph ),
	dot_write( DotF, Graph ).
