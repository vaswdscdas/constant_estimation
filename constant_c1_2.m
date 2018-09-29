%-----------------------------------------------------------------------
% Calcualte the P1 interpolation constant over triangle T, which 
% has the edges as e1,e2, e3.
% 
% $\Pi u$ is a linear function satisfying
% \int_{e_i} \Pi u -u ds = 0, i=1,2,3
%
% The interpolation error estimation with constant C:
% || u - \PI u || <= C^0 |u|_2 for u in H^2(T)
% | u - \PI u |_1 <= C^1 |u|_2 for u in H^2(T)
%
% Notation
% C^0: C_1^{(2,0)}(T)
% C^1: C_1^{(2,1)}(T)
%
%-----------------------------------------------------------------------
%
% Parameters of functions:
%
% mpath: path of mesh files.
% norm_idx: 0 or 1. (C^0： norm_idx=0, C^1: norm_idx=1）
%
% Xuefeng LIU 2014/04/28/

% Check before publishing [in process]
% Started on 2018/05/08
% Interval version: 2018/09/29


function [c_value]=constant_c1_2(mpath,norm_idx)

    tri = load([mpath,'tri_n.dat']); node = load([mpath,'node.dat']);
    edge = load([mpath,'edge.dat']);
    domain = load([mpath,'domain.dat']);

    %The precision of node points, which will be used in rigourous computing.
    %This is only used for rigorous computing.    
    node_precision = 1E-14;
    node = I_set_interval_precision(node, node_precision);
    
    display(sprintf('node is interval?  %d \n', isintval(node)));
    
    h = get_max_edge_length(edge,node);    
    display(sprintf('mesh size: %s \n', I_sup(h)));

    
    node_num = size(node,1);

    [B,H, null_dof_idx ] = sub_m_constraint_p12(tri,node,edge,domain);
    [A0,A1,A2] = create_matrix_morley(tri,node,edge);

    n=size(B,2);

     P=sparse(eye(n,n)) - H*B;
     NV=P'*A2*P;
     if norm_idx == 0
        DV = P'*A0*P;
        Ch = (I_intval(0.074)*h)^2;
     else
        DV = P'*A1*P;
        Ch = (I_intval(0.1888)*h);
     end
     
    ind=1:n;
    ind(null_dof_idx)=[];     

    lambda = est_min_eig(A2(ind,ind),NB(ind,ind));
    
    %~ lambda = est_min_eig(NV,DV);
    lambda_lower = lambda/(1+(Ch)^2*lambda);
    
    c_value = 1./sqrt(lambda_lower);

end

function [OutM1, H, null_dof_idx] = sub_m_constraint(node,edge,domain)

    e_num = size(edge,1); n_num = size(node,1);
    num = n_num + e_num;
    OutM1 = zeros(3,num); 
    OutM2 = zeros(3,num); 
    null_dof_idx = [0,0,0];

    domain_edges = domain( [2,3,1],:) - domain( [1,2,3],:);

    for domain_edge_ind=1:3
        for k=1:e_num
            local_edge = edge(k,:);
            edge_nodes = node(local_edge,:);
            edge_vec = edge_nodes(2,:) - edge_nodes(1,:);
            
            if is_edge_on_bd_edge_of_tri_domain(edge_nodes, domain, domain_edge_ind)
                g_index = k;
                OutM1( domain_edge_ind, g_index ) = -1;
                if null_dof_idx(domain_edge_ind) == 0
                    OutM2( domain_edge_ind, g_index ) = 1; 
                    null_dof_idx(domain_edge_ind) = g_index;
                end
            end           
            
        end
    end    
    H=OutM2';
    
end
