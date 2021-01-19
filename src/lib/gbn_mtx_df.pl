
gbn_mtx_df_defaults( check_names(false) ).

/** gbn_mtx_df( -Mtx, +Df, +Opts ).
mtx_df( +Mtx, +Df ).
 
Mtx is the csv representation of R data frame Df.

Opts
  * check_names(ChNames=true)
     in R the default is also _true_
 
==
?- gbn_mtx_df( [row(a,b,c),row(1,2,3),row(4,5,6)], df1 ), 
   <- df1,
   csv_df( Csv1, df1 ).
 
$a
[1] 1 4
 
$b
[1] 2 5
  
$c
[1] 3 6

Csv1 = [row(a, b, c), row(1, 2, 3), row(4, 5, 6)].
  
==

This seems to get stuck for very large matrices, (>130,000).
See implementation in r_sqlite_load.pl .
 
@author nicos angelopoulos
@version  0.1 2014/6/19
@version  0.2 2015/12/14, added options (check_names/1)
*/
gbn_mtx_df( Csv, Df ) :-
	gbn_mtx_df( Csv, Df, [] ).

gbn_mtx_df( Csv, Df, _Args ) :-
	var( Csv ),
	!,
	Rlist <- as.list(Df),
	findall( List, (member(Head=Tail,Rlist),List=[Head|Tail]), Lists),
	mtx_lists( Csv, Lists ).

gbn_mtx_df( Csv, Df, Args ) :-
	options_append( gbn_mtx_df, Args, Opts ),
	% Df <- 'data.frame()',
	mtx_lists( Csv, Lists ),
	findall( Head=Tail, ( member(List,Lists), List=[Head|Tail] ), Pairs ),
	/*
	findall( Head=Fail, ( member(List,Lists), List=[Head|Tail],
					  break_nth( 136560, Tail, Fail, Rest ),
					  Rest = [Off|_], write( Off ), nl
	                    ), Pairs ),
					*/
	Df <- Pairs, % this actually makes it to a list (if that matters)
	( options(check_names(true),Opts) -> Check = 'T'; Check = 'F' ),
	Df <- 'data.frame'(Df,check.names=Check).

% This ???Fixes a bug in one of: Real, R, Swi  (see ~/ac/14mg/hmrn/analysis
%
gbn_mtx_df_1( Csv, Df ) :-
	% Df <- 'data.frame()',
	mtx_lists( Csv, Lists ),
	findall( Head=Tail, ( member(List,Lists),List=[Head|Tail], (maplist(number,Tail)->true; throw(erto(Head))) ), Pairs ),
	/*
	findall( Head=Fail, ( member(List,Lists), List=[Head|Tail],
					  break_nth( 136560, Tail, Fail, Rest ),
					  Rest = [Off|_], write( Off ), nl
	                    ), Pairs ),
					*/
	Df <- Pairs, % this actually makes it to a list (if that matters)
	Df <- 'data.frame'(Df).
