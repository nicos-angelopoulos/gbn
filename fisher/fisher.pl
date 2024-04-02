
:- lib(real).
:- debuc(fish).

/** fisher.

	Make definitive test cases for Fisher colours

	From the tests below i concluded that 

	ft.pval < 0.05 is significance assocation

	ft.estimate < 1 is mutual exclusivity 

	ft.estimate > 1 is co-occurance 

*/

fisher :-
	X = c(1,1,1,1,0,0,0,0),
	Y = c(0,0,0,0,1,1,1,1),
	Z = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),
	W = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1),
	T = c(1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0),
	A = c(1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0),
	repo( f1, X, Y ),
	repo( f2, X, X ),
	repo( f3, Z, Z ),
	repo( f4, Z, W ),
	repo( f5, Z, T ),
	repo( f6, Z, A ).

repo( Rvar, V1, V2 ) :-
	Rvar <- fisher.test(V1,V2),
	debuc( fish, 'Rvar: ~w', Rvar ),
	<- Rvar,
	Pval <- Rvar$p.val,
	Estim<- Rvar$estimate,
	debuc(fish, 'pval:~w, estim:~w', [Pval,Estim]).
