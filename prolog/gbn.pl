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
:- use_module(library(lists)).
:- use_module(library(filesex)).
:- use_module(library(process)).

user:file_search_path( cancer, pack('gbn/run/gbns_in_cancer') ).
:- multifile(user:display_var_as/2).

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
% :- lib(gbn_module/0).
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

/** <module> Running and managing runs of gobnilp.


Run a GOBNILP, Bayesian networks (BNs) learning task and/or post-run routines on learning output.

Post processing includes the generation of Fisher based visualisation and family heatmaps.<br>
Assumes gobnlip is in your path. Currently only tested on *nix based systems.

---++ Installation

To install
==
?- pack_install(gbn).
==

to load
==
?- use_module(library(gbn)).
==

---++ Predicates

Main predicate for both running the BN experiment and post-processing is:
  * gbn/1

---++ Examples

==
?- gbn.
true.

?- ls.
% asia-24.07.13/   
true.
==

==
?- gbn(debug(true)).
% Turning debugging on for predicate handle: gbn(gbn)
% Options: [$restore(gbn,debug,false),copy(false),data(pack(gbn/data/asia.dat)),display_dot(svg),odir(_33202),std_output(std_file)]
% Output directory: 'asia-24.07.13'
% Settings on: asia.set
true.
==

A more complex example

==
?- [cancer(aml)].

?- absolute_file_name( pack('gbn/run/gbns_in_cancer'), Abs ),
|    ls( Abs ).
% aml.pl    coa.pl    data/     gbm.pl    mpn.pl    mye.pl    plots/    
Abs = '/home/nicos/.local/share/swi-prolog/pack/gbn/run/gbns_in_cancer'.

?- aml.
% Starting: aml
% Starting: fisher_nets
% Starting: fam_hmaps
% Starting: gates_nets
% Starting: svg_legend
% Finished: aml
true.
==

---++ Info
@author nicos angelopoulos
@version  0.0.1 2014/4/8
@version  0.1.0 2021/1/19
@version  0.2.0 2021/1/23
@see https://doi.org/10.1038/s42003-022-03243-w
@see https://stoics.org.uk/~nicos/sware/gbn
@see https://www.cs.york.ac.uk/aig/sw/gobnilp/
@see gbn/1
@see gbn_version/2

*/

/** gbn_module.

Documentation predicate to give anchor to the overall module documentation.

*/
gbn_module.

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
