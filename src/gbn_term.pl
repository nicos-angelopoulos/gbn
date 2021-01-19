
% :- lib(stoics_lib:io_lines/2).
% :- lib(stoics_lib:break_on_list/4).

/** gbn_term( +BnF, -Bn, -Adorned ).

Parse a .bn gobnilp file into a Bn term of the form [Ch-Pas|_] and 
Adorned BN that also includes the overall score and a score for each Ch-Pas entry
[Ch-Pas-Score|_]-BnScore.

==
?- absolute_file_name( pack(gbn/examples/crc.bn), Crc ), 
   gbn_term( Crc, Bn, Adorned ).

Crc = '/home/na11/lib/swipl/pack/gnb/examples/crc.bn'.
Bn = ['ACVR2A'-['ARID1A', 'FBXW7', 'PTEN'], 'AMER1'-[], 'APC'-['ACVR2A', 'AMER1'], 'ARID1A'-[], 'ATM'-['ACVR2A'], 'BRAF'-['APC'], 'FBXW7'-[], 'KRAS'-[...], ... - ...|...],
Adorned = ['ACVR2A'-['ARID1A', 'FBXW7', 'PTEN']- -60.029495, 'AMER1'-[]- -79.527071, 'APC'-['ACVR2A', 'AMER1']- -131.057359, 'ARID1A'-[]- -84.299473, 'ATM'-['ACVR2A']- -43.772266, 'BRAF'-[...]- -110.700005, ... - ... - -86.616242, ... - ...|...]- -1567.186608.
==
*/

gbn_term( BnF, Bn, Abn-Score  ) :-
    io_lines( BnF, BnLines ),
    once( append(FamLines,[ScoreLine],BnLines) ),
    gbn_score_line( ScoreLine, Score ),
    maplist( gbn_line_family, FamLines, Bn, Abn ).

gbn_line_family( Line, Ch-Pas, Ch-Pas-Score ) :-
    break_on_list( Line, [0'<,0'-], ChCs, PasCs ),
    name( Ch, ChCs ),
    bn_comma_seperated_parents( PasCs, Pas, Score ).

bn_comma_seperated_parents( Codes, Pas, Score ) :-
    break_on_list( Codes, [0',], Left, Right ),
    !,
    name( Pa, Left ),
    Pas = [Pa|TPas],
    bn_comma_seperated_parents( Right, TPas, Score ).
bn_comma_seperated_parents( PrvCodes, [], Score ) :-
    ( PrvCodes = [0' |Codes] -> true; Codes = PrvCodes ),
    name( Score, Codes ).

gbn_score_line( Line, Score ) :-
    break_on_list( Line, [0' ,0'i,0's,0' ], _Left, Right ),
    number_codes( Score, Right ).
