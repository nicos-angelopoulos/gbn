/** gbn_version( -Vers, -Date ).

Get the current version and publication date.

==
?- gbn_version( Vers, Date ).
V = 0:2:0,
D = date(2021, 1, 23).

==

@author nicos angelopoulos
@version  0.0.1 2014/4/8
@version  0.1.0 2021/1/19
@version  0.2.0 2021/1/23
@see http://stoics.org.uk/~nicos/sware/gbn
@see http://www.cs.york.ac.uk/aig/sw/gobnilp/
@tbd control replacement of dot (and other specials?) in Bn var names via options

*/
gbn_version(0:2:0, date(2021,1,23)).

gbn_defaults( Defs ) :-
	Defs = [ % dir(Gob),
              copy(false),
		    data(pack('gbn/data/asia.dat')),
		    display_dot(svg),
	         std_output(std_file)
	       ].

gbn :-
	gbn( [] ).

/** gbn( Options ).

Options is a term or list of:

  * copy(Copy=false)
    if set to anything else than false the local file Copy is copied into Dir- multiple are allowed

  * dir(Dir)
    output directory for runs, it should not exist (default: *<Pbname>-<date>*). When a variable, the name is returned

  * display_dot(DispDot=svg)
    display dot filin format DispDot. false for not displaying.

  * data(Dfile) 
    The datafile on which gobnilp will be ran (*pack(gbn/data/asia_100.dat*)).

  * debug(Dbg)  
    can be =|false|= (also _true_=_gbn_ or _all_).

  * multiple(Setting,Vals)
    run multiple experiments where setting is given each Val-ue in turn

  * probname(Pbname)
    used as stem for output files. Defaults is the stem of Dfile (so for Dfile =../data/asia.dat the stem is _asia_).

  * std_ouput(StdO)   
    what to do with gobnilp output ?, (*std_file*,output,<file>). Where std_file is <Pbname>_std_output.txt

  * setting(Key,Val)
    set Gobnilp setting Key to value Val

@author nicos angelopoulos
@version  0.1 2014/4/8
@tbd  subdirs option (ie non flat run dirs)
@tbd  bootstraps
@tbd  hook on to the c-code  (see pack(gob))

*/
gbn( Args ) :-
	options_append( gbn, Args, All, [process(debug),debug(true),module(gbn)] ),
	debuc( gbn(gbn), 'Options: ~w', [All] ),
	selectchk( data(Dfile), All, NtAll ),
	absolute_file_name( Dfile, AbsDfile ),
	file_base_name( AbsDfile, DBname ),

	gbn_prob_name( NtAll, ProbName, DBname, NpAll ),

	gbn_out_dir( ProbName, Dir, All ),
	make_directory( Dir ),

	directory_file_path( Dir, DBname, DatTarget ),
	copy_file( AbsDfile, DatTarget ),
    findall( _, (member(copy(Copy),All),
                   ( Copy == false -> true; @cp(Copy,Dir) )
                ),
                    _ ),
	% rename_file( SetsF, Dir ),

	debuc( gbn(gbn), 'Output directory: ~p', [Dir] ),
	working_directory( Old, Dir ),
	gob_prob_name_settings( ProbName, DBname, SetsF, All ),
	gbn_multiples( ProbName, DBname, SetsF, NpAll ),

	working_directory( _, Old ).

gbn_singleton( ProbName, DBname, SetsF, All ) :-
	debuc( gbn(gbn), 'Settings on: ~w', SetsF ),
	atomic_list_concat( ['-g=',SetsF], Garg ),
	selectchk( std_output(Sout), All, _NsAll ),
	std_output_gobnilp( Sout, Garg, ProbName, DBname ),
	options( display_dot(DispDot), All ),
	gbn_singleton_display_dot( DispDot, SetsF ).

gbn_singleton_display_dot( false, _SetsF ) :-
	!.
% fixme: allow Ext that are not the same as their format
gbn_singleton_display_dot( Ext, SetsF ) :-
	os_ext( set, dot, SetsF, DotF ),
	os_ext( set, Ext, SetsF, ImgF ),
	atom_concat( '-T', Ext, TForm ),
	% gobnilp does not quote specials, so HLA-A trips the dot exec.
	% pupsh is my local script for fixing those files by replacing with "HLA-A"
	% maybe we can do such calls via options ... or hooks
	% ( once(which(pupsh,_)) ->
    dot_fix( DotF ),
    ( which(dot,_DotEx) ->
	    @ dot( TForm, -o, ImgF, DotF )
        ; 
        % fixme: create info channel ?
        gbn_message( exec_miss(dot,'cannot create net',DotF) )
    ).

gbn_multiples( ProbName, DBname, SetsF, All ) :-
	\+ memberchk( multiple(_,_), All ),
	!,
	gbn_singleton( ProbName, DBname, SetsF, All ).
gbn_multiples( ProbName, DBname, SetsF, All ) :-
	% gbn_singleton( ProbName, DBname, SetsF, All ),
	findall( t(Key,Ktkn,Vals), 
					( member(multiple(Key,Vals),All),
					  ( memberchk(setting_token(Key,Ktkn),All) ->
					  	true
						;
						downcase_atom(Key,Ley),
						sub_atom(Ley,0,1,_,Ktkn)
					  )
					), KKVss ),
	exclude( form_term(multiple/2), All, Opts ),
	once((gbn_multiples_loop(KKVss,ProbName,DBname,SetsF,'',Opts);true)).

gbn_multiples_loop( [], ProbName, DBname, SetsF, StemExt, All ) :-
	os_postfix( StemExt, SetsF, ToSetsF, sep(-) ),
	@ cp( -i, SetsF, ToSetsF ),
	gbn_settings_file_append( ToSetsF, All ),
	gbn_singleton( ProbName, DBname, ToSetsF, [display_dot(false)|All] ),
	% fixme: change svg, to DispDot option
	maplist( gbn_multiple_rename(SetsF,ToSetsF), [dot,bn] ),
	options( display_dot(DispDot), All ),
	gbn_singleton_display_dot( DispDot, ToSetsF ).

gbn_multiples_loop( [t(Key,Tkn,Vals)|T], ProbName, DBname, SetsF, StemExt, Opts ) :-
	member( Val, Vals ),
	atomic_list_concat( [StemExt,Tkn,Val], NxtExt ),
	All = [setting(Key,Val)|Opts],
	gbn_multiples_loop( T, ProbName, DBname, SetsF, NxtExt, All ),
	fail.

gbn_settings_file_append( SetsF, Opts ) :-
	open( SetsF, append, Out ),
	findall( _, (member(setting(Key,Val),Opts),gob_stream_setting(Out,Key,Val)), _ ),
	close( Out ).

form_term( Name/Arity, Term ) :-
	functor( Term, Name, Arity ).

gbn_multiple_rename( SetsF, ToSetsF, Ext ) :-
	os_ext( Old, Ext, SetsF, From ),
	os_ext( Old, Ext, ToSetsF, To ),
	@ mv( -i, From, To ).

gbn_out_dir( _ProbName, Dir, All ) :-
	memberchk( dir(Dir), All ),
	ground( Dir ),
	!.
gbn_out_dir( _ProbName, Dir, Opts ) :-
	memberchk( dir_prefix(Pfx), Opts ),
	!,
	os_unique( Pfx, Dir, [] ).
gbn_out_dir( ProbName, Dir, Opts ) :-
	os_unique( ProbName, Dir, [create(false)] ),
	( memberchk(dir(Dir),Opts) -> true % returns it
	                            ; true ).

std_output_gobnilp( output, Garg, _Prob, DBname ) :-
	@ gobnilp( Garg, DBname ).
std_output_gobnilp( std_file, Garg, Prob, DBname ) :-
	atom_concat( Prob, '_std_output', ProbStdOut ),
	file_name_extension( ProbStdOut, txt, PSOf ),
	!,
	file_ouput_gobnilp( PSOf, Garg, DBname ).

file_ouput_gobnilp( OutF, Garg, DBname ) :-
	Outs @@ gobnilp( Garg, DBname ),
	open( OutF, write, Out ),
	maplist( writenl(Out), Outs ),
	close( Out ).

writenl( Out, Line ) :-
	write( Out, Line ), nl( Out ).

% file_name_extension( ProbName, bn, BnF ),

gbn_prob_name( All, ProbName, _DBname, Red ) :-
	select( probname(ProbName), All, Red ),
	!.
gbn_prob_name( All, ProbName, DBname, All ) :-
	file_name_extension( ProbName, _Est, DBname ).
	
gob_prob_name_settings( ProbName, _DBname, SetsF, Opts ) :-
	file_name_extension( ProbName, set, SetsF ),
	open( SetsF, write, Out ),
	gob_stream_setting( Out, 'nbns', 1 ),
	gob_stream_setting( Out, 'outputfile/solution', ProbName, bn ),
	% gob_stream_setting( Out, 'outputfile/solution', '"ran-asia/asia.bn"' ),
	gob_stream_setting( Out, 'outputfile/dot', ProbName, dot ),
	atom_concat( ProbName, '_scnti', Scnti ),
	gob_stream_setting( Out, 'outputfile/scoreandtime', Scnti, txt ),
	findall( _, (member(setting(Key,Val),Opts),gob_stream_setting(Out,Key,Val)), _ ),
	close( Out ).

gob_stream_setting( Out, Key, Value ) :-
	write( Out, 'gobnilp/' ), 
	write( Out, Key ), 
	write( Out, ' = ' ),
	( number(Value) -> write(Out,Value); 
				    write(Out,'"' ),
				    write(Out,Value),
				    write(Out,'"' ) ),
	nl( Out ).
gob_stream_setting( Out, Key, Value, Ext ) :-
	write( Out, 'gobnilp/' ),
	write( Out, Key ),
	write( Out, ' = "' ),
	write( Out, Value ),
	write( Out, '.' ),
	write( Out, Ext ),
	write( Out, '"' ),
	nl( Out ).
