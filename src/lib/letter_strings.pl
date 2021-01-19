/** letter_strings( +Start, -N, -Letts ).

Generate N letter strings, starting from Start.
Start is polymorphic: string, code (integer) or atom.

==
?- letter_strings( a, 3, Letts ).
Letts = ["a", "b", "c"].

?- letter_strings( "C", 3, Letts ).
Letts = ["C", "D", "E"].

==
@author nicos angelopoulos
@version  0.1 2017/2/15
@tbd      check we do not over-run

*/
letter_strings( Start, N, Letts ) :-
    atom( Start ), 
    !,
    atom_codes( Start, [Code] ),
    letter_strings_1( N, Code, Letts ).
letter_strings( Start, N, Letts ) :-
    integer( Start ), 
    !,
    letter_strings_1( N, Start, Letts ).
letter_strings( Start, N, Letts ) :-
    string( Start ),
    !,
    string_codes( Start, [Code] ),
    letter_strings_1( N, Code, Letts ).
letter_strings( [Code], N, Letts ) :-
    letter_strings_1( N, Code, Letts ).

letter_strings_1( 0, _Code, [] ) :- !.
letter_strings_1( N, Code, [L|Ls] ) :-
    string_codes( L, [Code] ),
    M is N - 1,
    Node is Code + 1,
    letter_strings_1( M, Node, Ls ).
