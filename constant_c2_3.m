%-----------------------------------------------------------------------
% Calcualte the Fujino-Morley interpolation constant over triangle T, which 
% has the vertices as p1,p2,p3, and edges as e1,e2,e3.
% 
% $\Pi u$ is a quadratic polynomial satisfying
% (\Pi u) (p_i) = u(p_i), i=1,2,3
% \int_{e_i} \partial (\PI u - u)/\partial n ds = 0, i=1,2,3
%
% The interpolation error estimation with constant C:
% || u - \PI u || <= C^0 |u|_2 for u in H^2(T)
% | u - \PI u |_1 <= C^1 |u|_2 for u in H^2(T)
%
% In case norm_idx = 0:  c = C^0
% In case norm_idx = 1:  c = C^1
%-----------------------------------------------------------------------
%
% Xuefeng LIU 2012/10/05
% Check before publishing [in process]
% Started on 2018/05/08
% Interval version: 2018/09/26

function [c_value]=constant_c2_3(mpath,norm_idx)

    tri = load([mpath,'tri_n.dat']); node = load([mpath,'node.dat']);
    edge = load([mpath,'edge.dat']);
    edge = sort(edge,2);
    domain = load([mpath,'domain.dat']);

    ind1=find( ([abs(node(:,1)-domain(1,1) ), abs(node(:,2)-domain(1,2))]*[1,1]') < 1E-6 );
    ind2=find( ([abs(node(:,1)-domain(2,1)), abs(node(:,2)-domain(2,2))]*[1,1]') < 1E-6 );
    ind3=find( ([abs(node(:,1)-domain(3,1)), abs(node(:,2)-domain(3,2))]*[1,1]') <1E-6 );

    %The precision of node points, which will be used in rigourous computing.
    %This is only used for rigorous computing.    
    node_precision = 1E-14;
    node = I_set_interval_precision(node, node_precision);
    
    display(sprintf('Is node in interval mode?  %d \n', isintval(node)));
    
    h = get_max_edge_length(edge,node);    
    display(sprintf('Mesh size: %s \n', I_sup(h)));


    [C,H, null_dof_idx]=sub_m_constraint(node,edge,domain);
    
    [A0,A1,A2]=create_matrix_morley(tri,node,edge);

    % % Method 1
    % 
    % NA=[A2, B'; B, zeros(3,3)];
    % NB=[A0, 0*B'; 0*B, zeros(3,3)];
    % 
    % n=size(NA,1);
    % ind=1:n;
    % ind([ind1,ind2,ind3])=[];
    % 
    % DDV = NA(ind,ind);
    % DV = NB(ind,ind);
    % [v,d] = eigs(DDV, DV, 1, 'sm');
    % c_value = 1./sqrt(d);


    % Method 2
    n=size(A0,1);
    P=sparse(eye(n,n)) - H*C;

    if norm_idx == 0
        NB = P'*A0*P;
        Ch = (I_intval(0.0736)*h)^2;
    else
        NB = P'*A1*P;
        Ch = (I_intval(0.1887)*h);
    end

    ind=1:n;
    null_idx = [ind1,ind2,ind3, null_dof_idx]
    ind(null_idx)=[];     

    lambda = est_min_eig(A2(ind,ind),NB(ind,ind));
    lambda_lower = lambda/(1 + Ch^2*lambda);
    
    c_value =1./sqrt(lambda_lower);

    return

    figure(2)
    hold off
    node_num = size(node,1);
    vec=zeros(node_num,1); vec(ind)=v(:,1);
    trisurf(tri,node(:,1), node(:,2),vec(1:node_num));

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
                if edge_vec*domain_edges(domain_edge_ind,:)'>0
                    edge_sign = 1;
                else
                    edge_sign = -1;
                end
                g_index = n_num + k;
                OutM1( domain_edge_ind, g_index ) = edge_sign;
                if null_dof_idx(domain_edge_ind) == 0
                    OutM2( domain_edge_ind, g_index ) = edge_sign; 
                    null_dof_idx(domain_edge_ind) = g_index;
                end
            end           
            
        end
    end    
    H=OutM2';
    
end
