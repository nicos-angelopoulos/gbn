/** dot_file_bn_fams( +DotF, -Fams ).

Create Ch-Pa Fam-ily representation of a BN defined in DotF.

This used to be pinched from .bn files via gbn_term/3, but recent GOBNILPs do not produce these.

@author nicos angelopoulos
@version  0.1 2026/05/26

*/

dot_file_bn_fams( DotF, Fams ) :-
     dot_read( DotF, dot(_,_,_,Graph) ),
     % GOBNILP seems to have all nodes first anyway, 
     % however here we look for them in edges as well
     findall( SiNode, (member(node(QSiNode,_),Graph),dot_file_deq(QSiNode,SiNode)), SiNodes ),   
               % ^ fixme also look for node(Node) ? 
     findall( [FromNode,ToNode], ( member(edge([QFromNode,'->',QToNode],_),Graph),
                                    dot_file_deq(QFromNode,FromNode),
                                    dot_file_deq(QToNode,ToNode)
                                 ), NestNodes ),
     flatten( [SiNodes,NestNodes], RepNodes ),
     list_to_ord_set( RepNodes, Nodes ),
     findall( Node-Pas, (member(Node,Nodes), findall(Pa,member(edge([Pa,'->',Node],_),Graph),Pas)), Fams ).

dot_file_deq( quote(Node), Node ) :- !.
dot_file_deq( Node, Node ).

