%-----------------------------------------------------------------------
% Calcualte the Crouzeix-Raviart interpolation constant  over triangle T, which 
% has the edges as e1,e2, e3.
% 
% $\Pi u$ is a linear function satisfying
% \int_{e_i} \Pi u -u ds = 0, i=1,2,3
%
% The interpolation error estimation with constant C:
% || u - \PI u || <= C | u - \PI u |_1 for u in H^1(T)
%
%-----------------------------------------------------------------------
%
% Parameters of functions:
%
% mpath: path of mesh files.
%
% Xuefeng LIU 2014/04/28/
%
% Check before publishing [in process]
% Started on 2018/05/08
% Interval version: 2018/09/26

function [c_value]=constant_crouzeix_raviart(mpath)

    tri = load([mpath,'tri_n.dat']); node = load([mpath,'node.dat']);
    edge = load([mpath,'edge.dat']);
    domain = load([mpath,'domain.dat']);


    t_num = size(tri,1);  e_num = size(edge,1);
    tri_by_edge = zeros(t_num,3); %Represent a triangle by 3 edges.
    edge_value_index = sort(edge,2)*[e_num;1]; %Create values for all edges, which helps to find an edge by its two nodes.

    for k=1:t_num  %create the table of element-edges relation.
        for l = 1:3  %the l-th edge.
            node_ind =[ mod(l,3)+1, mod(l+1,3)+1];
            value = sort(tri(k, node_ind) )*[e_num;1];
            [r,ind] = ismember( value, edge_value_index);
            tri_by_edge(k,l) = ind; %The direction of an edge is not counted, which is needed for other elements, e.g., Fujino-Morley FEM
        end
    end

    %The precision of node points, which will be used in rigourous computing.
    %This is only used for rigorous computing.    
    node_precision = 1E-14;
    node = I_set_interval_precision(node, node_precision);
    
    display(sprintf('node is interval?  %d \n', isintval(node)));
    
    h = get_max_edge_length(edge,node);    
    display(sprintf('mesh size: %s \n', I_sup(h)));

    [B,H, null_dof_idx ] = sub_m_constraint(node,edge,domain);
    [A0,A1] = create_matrix_crouzeix_raviart(tri,edge,node,tri_by_edge);
   
    n=e_num;

    P=sparse(eye(n,n)) - H*B;
    NA = P'*A1*P;
    NB = P'*A0*P;
    Ch = I_intval(0.1893)*h;

    ind=1:n;
    ind(null_dof_idx)=[];     

    lambda = est_min_eig(NA(ind,ind),NB(ind,ind));
    lambda_lower = lambda/(1+(Ch)^2*lambda);
    
    c_value = 1./sqrt(lambda_lower);

end

function [OutM1, H, null_dof_idx] = sub_m_constraint(node,edge,domain)

    e_num = size(edge,1); n_num = size(node,1);
    num = e_num;
    OutM1 = zeros(3,num); 
    OutM2 = zeros(3,num); 
    null_dof_idx = [0,0,0];

    for domain_edge_ind=1:3
        for k=1:e_num
            local_edge = edge(k,:);
            edge_nodes = node(local_edge,:);
            
            if is_edge_on_bd_edge_of_tri_domain(edge_nodes, domain, domain_edge_ind)
                g_index = k;
                OutM1( domain_edge_ind, g_index ) = -1;
                if null_dof_idx(domain_edge_ind) == 0
                    OutM2( domain_edge_ind, g_index ) = -1; 
                    null_dof_idx(domain_edge_ind) = g_index;
                end
            end           
            
        end
    end    
    H=OutM2';
    
end

