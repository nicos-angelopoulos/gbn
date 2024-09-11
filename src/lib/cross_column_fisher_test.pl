:- lib(real).

/** cross_column_fisher_test( +I, +J, +Cols, +Data, +SimPv, +RInters, +ROdds ).

Add pval and estimate (oods.ratio) of fisher.test on two columns
in Data to the R metrices RInters (interactions) and ROdds, respectively.

SimPv is an Real passable boolean ('TRUE'/'FALSE'), that is passed to fisher.test() as simulate.p.value. 

The pvalue of the test shows if there is an interaction and the odds 
show the type of intaction ( <1 is a mutually exclusive interaction,
and >1 is a co-occurance interaction ).

This is meant, and only-tested-on, for binary vectors.

The old version is commented out. It is likely it was seriously wrong.

@version 0:2  2016/12/5
@version 0:3  2024/4/2,  added SimPv

*/
cross_column_fisher_test( I, J, Cols, Data, SimPv, Inters, Odds ) :-
	I > Cols,
	NxJ is J + 1,
	!,
	cross_column_fisher_test_nest( NxJ, Cols, Data, SimPv, Inters, Odds ).
cross_column_fisher_test( I, J, Cols, Data, SimPv, Inters, Odds ) :-
	% debug( gbn(fisher), 'I: ~d, J: ~d', [I,J] ),
	f <- try(fisher.test(Data[*,I], Data[*,J],simulate.p.value=SimPv), silent='TRUE'),
	% fisher_log10_odds( f, Log10Val, OddVal ),
	fisher_pval_odds( f, Pval, OddVal ),
	% Inters[I,J] <- Log10Val,
	Inters[I,J] <- Pval,
	IsInf <- is.infinite(OddVal),
	( IsInf == true ->
          Left <- 'f$conf.int[1]',
          ( Left >= 1 -> % fixme: not sure 1 is best... exclusive seem to be [0.0000,Inf] ; co-occur seem [8k,Inf]
		     Odds[I,J]   <- 1/0
               ;
		     Odds[I,J]   <- -1/0   % fixme, on prompt -1.0Inf does not work on SWI prompt
          )
		% Odds[I,J]   <- 'NA'
		;
		Odds[I,J]   <- OddVal
	),
	NxI is I + 1,
	cross_column_fisher_test( NxI, J, Cols, Data, SimPv, Inters, Odds ).
	
cross_column_fisher_test_nest( J, Cols, _Data, _SimPv, _Inters, _Odds ) :-
	J > Cols,
	!.
cross_column_fisher_test_nest( J, Cols, Data, SimPv, Inters, Odds ) :-
	cross_column_fisher_test( 1, J, Cols, Data, SimPv, Inters, Odds ).

fisher_pval_odds( Fv, Pval, OddVal ) :-
	Error <- class(Fv),
	( Error == 'try-error' ; (Est <- Fv$estimate,  Est == []) ), % double check the RHS of ; (own addition)
	!,
	OddVal = 'NA',  % check
	Pval is 1.
fisher_pval_odds( Fv, Pval, OddVal ) :-
	OddVal = Fv$estimate,
	Pval <- Fv$p.val.

/* I think (16.12.05 this is wrong as it logs the pval ?!?
   if there is any need to log things that would be 

fisher_log10_odds( Fvar, Log10Val, OddVal ) :-
	Error <- class(Fvar),
	( Error == 'try-error' ; (Est <- Fvar$estimate,  Est == []) ), % double check the RHS of ; (own addition)
	!,
	Log10Val is 0,
	OddVal = 'NA'.  % check
fisher_log10_odds( Fvar, Log10Val, OddVal ) :-
	Est <- Fvar$estimate,
	Est > 1,
	!,
	PvalPrv <- Fvar$p.val,
	% debug( _, 'Pval: ~w', PvalPrv ),
	( PvalPrv =:= 0 -> Pval is 10e-10; Pval = PvalPrv ),
	Log10Val is - log10(Pval),
	OddVal = Fvar$estimate.
	% OddVal <- Fvar$estimate.
fisher_log10_odds( Fvar, Log10Val, OddVal ) :-
	Pval <- Fvar$p.val,
	Log10Val is log10(Pval),
	% OddVal <- Fvar$estimate.
	OddVal = Fvar$estimate.
	*/
