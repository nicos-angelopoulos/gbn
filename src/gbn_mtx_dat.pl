
:- lib(mtx).
:- lib(suggests(os_lib)).

/** gbn_mtx_dat( +CsvF, ?DatF ).
	
	Convert csv input file Csv to a space separated value file
	which contains as a second line the counts of distinct values for each column.
	If DatF is a variable it is instantiatiated to file indicator which is produced by  replacing
	the extension in CsvF with '.dat' . 
	These .dat files are input data files to gobnilp.

	Requires pack(mtx) and pack(os) to work. pack(os) is only required if DatF is a variable.
	Loading does not complain if this lib does not load properly (this is to avoid annoying message
	for user that do not require this modality of the predicate).

@author nicos angelopoulos
@version  0.1 2016/1/27
@see mtx/2       CsvF can be any recognised matrix input to mtx/2
@see suggests/1  (pack(requires))

*/
gbn_mtx_dat( CsvF, DatF ) :-
	mtx( CsvF, Mtx ),
	mtx_columns_values( Mtx, Csets, [has_header(true),values_as(set)] ),
	maplist( length, Csets, Lengths ),
	( var(DatF) -> os_ext(_,dat,CsvF,DatF); true ),
	mtx_header_body( Mtx, Hdr, Body ),
	PopR =.. [row|Lengths],
	csv_write_file( DatF, [Hdr,PopR|Body], [separator(0' )] ).
