:- module( gbn, [   gbn/0,  gbn/1,  
                    gbn_version/2, gbn_module/0,
                    gbn_family_gates/5,
                    gbn_fam_hmaps/1,
                    gbn_fisher_net/3, gbn_fisher_net/4, 
                    gbn_fisher_nets/0,gbn_fisher_nets/1,
                    gbn_mtx_dat/2,gbn_mtx_paired/4,gbn_mtx_subs/3,
                    gbn_mtx_mins/3,
                    gbn_term/3,
                    gbn_prefix_directed_constraints/4,
                    gbn_gates_nets/0, gbn_gates_nets/1,
                    gbn_gates_net/1,gbn_gates_net/2,
                    gbn_svg_legend/1
                    ]  ).

:- use_module(library(csv)).
:- use_module(library(pio)).
:- use_module(library(apply)).
:- use_module(library(filesex)).
:- use_module(library(process)).

:- use_module(library(lib)).
:- lib(source(gbn), homonyms(true)).

:- lib(mtx).
:- lib(real).
:- lib(os_lib).
:- lib(by_unix).
:- lib(disp_bn).
:- lib(options).
:- lib(debug_call).
:- lib(stoics_lib).
:- lib(suggests(svg)).

:- lib(r(ggplot2)).

:- lib(gbn/1).
:- lib(gbn_module/0).
:- lib(gbn_term/3).
:- lib(gbn_fisher_net/3).
:- lib(gbn_fisher_net/4).
:- lib(gbn_fisher_nets/0).
:- lib(gbn_fisher_nets/1).
:- lib(gbn_mtx_dat/2).
:- lib(gbn_mtx_subs/3).
:- lib(gbn_fam_hmaps/1).
:- lib(gbn_mtx_paired/4).
:- lib(gbn_family_gates/5).
:- lib(gbn_prefix_directed_constraints/4).
:- lib(gbn_gates_net/1).
:- lib(gbn_gates_net/2).
:- lib(gbn_gates_nets/0).
:- lib(gbn_gates_nets/1).
:- lib(gbn_svg_legend/1).
:- lib(gbn_mtx_mins/3).
:- lib(gbn_res_dir_dat_file/2).
:- lib(dot_fix/1).
:- lib(gg_muts_by_pnms/3).

:- lib(end(gbn)).

user:file_search_path( cancer, pack('gbn/run/gbns_in_cancer') ).


                 /*******************************
                 *            MESSAGES          *
                 *******************************/
% These print messages that are always on.
% Different colour to debugging is used by the system (when colour in terminal is enabled).
%
gbn_message( Mess ) :-
    print_message( informational, gbn(Mess) ).
    
:- multifile prolog:message//1.

prolog:message(gbn(Message)) -->
    message(Message).

:- discontiguous
    message//1.

message( exec_miss(Exec,Facil,Target) ) -->
    ['The ~w executable is not installed, ~w for ~w'-[Exec,Facil,Target] ].
message( non_unique_dat([],Dir) ) -->
    ['Directory ~w does not contain the needed .dat file.'-[Dir] ].
message( non_unique_dat([_A,_B|_T]), Dir ) -->
    ['Directory ~w contains more than one .dat files.'-[Dir] ].
message( cow_plot_lbl(Tkn,Use) ) -->
    ['Wrong cow plot label: ~w, using: ~w.'-[Tkn,Use] ].
