
gbn_gates_nets_defaults( dir('.') ).

/** gbn_gates_nets.
    gbn_gates_nets( +Opts ).

Run gbn_gates_net/1 on all dot files in Dir (default is '.').

Opts 
  * dir(Dir='.')
     set working directory

@author nicos angelopoulos
@version  0.1 2017/09/22

*/
gbn_gates_nets :-
    gbn_dot_gates_dir( [] ).

gbn_gates_nets( Args ) :-
    options_append( gbn_gates_nets, Args, Opts ),
    options( dir(Dir), Opts ),
    os_sel( os_files, ext(dot), DotFs, dir(Dir) ),
    include( os_postfix(fclr), DotFs, GatFs ),
    maplist( os_gates_dir_net(Dir), GatFs ).

os_gates_dir_net( Dir, DotF ) :-
    os_path( Dir, DotF, DotP ),
    gbn_gates_net( DotP ).
