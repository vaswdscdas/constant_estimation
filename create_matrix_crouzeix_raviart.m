function [A0, A1] = create_matrix_crouzeix_raviart(tri, edge, node, tri_by_edge)
%% Function to calcuate the inner product of P1 bases of Crouzeix-Ravart FEM
%% Xuefeng LIU (2014/06/16)
%% An interval version can be found in <Workspace>/FEM_Code_LIB/create_matrix_crouzeix_raviart.m

%     tri = load('tri.dat'); node = load('node.dat');
%     tri = load([mpath,'tri_n.dat']); node = load([mpath,'node.dat']);    

    t_num = size(tri,1); n_num = size(node,1);  e_num = size(edge,1);  num = e_num;

    A0 = sparse( num, num );  
    A1 = sparse( num, num );

    for k=1:t_num
        t=tri(k,:);
        edges = node( t([3,1,2]), : )  - node( t([2,3,1]), : ); %3 by 2
        edge_index = tri_by_edge(k,:);        
        S =  abs(0.5* edges(1,:)*[edges(2,2); -edges(2,1)]);
        e_e = edges*edges';
        A0( edge_index, edge_index ) = A0( edge_index, edge_index ) + S*eye(3,3)/3.0;
        A1( edge_index, edge_index ) = A1( edge_index, edge_index ) + e_e/S;
    end
end
