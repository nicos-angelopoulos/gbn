/** gbn_bn_neighs( +Bn, ?Node, -Pa, -Ch ).

True if Pa is the set of parents for Node in Bn, and Ch
is the set of children.

==
?- requires( gbn:gbn_bn_neighs/4 ).
==
*/

gbn_bn_neighs( Bn, Node, Pas, Chs ) :-
	( var(Node) ->
		member( Node-Pas, Bn )
		;
		memberchk( Node-Pas, Bn )
	),
	% findall( )
	findall( Ch, (member(Ch-ChPas,Bn),memberchk(Node,ChPas)), ChsL ),
	sort( ChsL, Chs ).
