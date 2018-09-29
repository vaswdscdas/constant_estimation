%-----------------------------------------------------------------------
% Calcualte the Fujino-Morley interpolation constant over triangle T, which 
% has the vertices as p1,p2,p3, and edges as e1,e2,e3.
% 
% $\Pi u$ is a quadratic polynomial satisfying
% (\Pi u) (p_i) = u(p_i), i=1,2,3
% \int_{e_i} \partial (\PI u - u)/\partial n ds = 0, i=1,2,3
%
% The interpolation error estimation with constant C:
% || u - \PI u || <= C^0 |u - \PI u |_2 for u in H^2(T)
% | u - \PI u |_1 <= C^1 |u - \PI u |_2 for u in H^2(T)
%
% In case norm_idx = 0:  c = C^0
% In case norm_idx = 1:  c = C^1
%-----------------------------------------------------------------------
%
% Xuefeng LIU started on 2012/10/05
%
% Check before publishing [in process]
% Started on 2018/05/08
% Interval version: 2018/09/26

function [c_value]=constant_fujino_morley(mpath,norm_idx)

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
    
    display(sprintf('Is node in interval node? %d \n', isintval(node)));
    
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

    NA = P'*A2*P;

    if norm_idx == 0
        NB = P'*A0*P;
        Ch = (I_intval(0.0736)*h)^2;
    else
        NB = P'*A1*P;
        Ch = (I_intval(0.1887)*h);
    end

    ind=1:n;
    null_idx = [ind1,ind2,ind3, null_dof_idx];
    ind(null_idx)=[];     

    lambda = est_min_eig(NA(ind,ind),NB(ind,ind));
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


%~ function [OutM1, H] = sub_m_constraint(tri,node,edge,domain)

%~ %     domain = load([mpath,'domain_tri.dat');
    %~ e_num = size(edge,1); t_num = size(tri,1); n_num = size(node,1);
    %~ tri_by_edge = zeros(t_num,3);
    %~ edge_value_index = edge*[n_num,1]';
    %~ num = n_num + e_num;
    %~ OutM1 = sparse( 3,num); 
    %~ OutM2 = sparse( 3,num); 
    %~ for k=1:t_num  %create the table of element-edges relation.
        %~ for l = 1:3  %the l-th edge.
            %~ node_ind =[ mod(l,3)+1, mod(l+1,3)+1];
            %~ local_value = sort( tri(k, node_ind) )*[n_num,1]';
            %~ ind = find( edge_value_index - local_value == 0);
            %~ tri_by_edge(k,l) = - ind * ( 2*( norm( tri(k,node_ind)-edge(ind,:) ) >0 ) - 1);
        %~ end
    %~ end
    
    %~ domain_edges = domain( [3,1,2],:) - domain( [2,3,1],:);
    %~ vec1vec2= [domain_edges(2,:)', -domain_edges(1,:)'];    
    %~ st_coord = vec1vec2\([node(:,1) - domain(3,1), node(:,2)-domain(3,2)]');
    %~ st_coord = st_coord';
    %~ % What we need to do?
    %~ % 1) three speical functions under area coordinates.
    %~ % 2) the value of functions on each nodes of local element.
    %~ % 3) the coefficent of functions at each local freedom.

    %~ d_S =  0.5* domain_edges(1,:)*[domain_edges(2,2); -domain_edges(2,1)];        
    %~ d_e_e = domain_edges*domain_edges';
    %~ d_C =  [-1,-1,1;1,-1,-1;-1,1,-1]*diag(diag(d_e_e))/4.0/d_S ;
    %~ coeff_Function = inv(d_C); % coefficient in each row

    
    %~ for k=1:t_num
        %~ t=tri(k,:);
        %~ edge_index = tri_by_edge(k,:); %this can be negative.
        %~ edge_index_sign = [1,1,1, (2*( edge_index>0 ) -1 )];
        %~ local_edges = node( t([3,1,2]), : )  - node( t([2,3,1]), : );
        %~ S =  0.5* local_edges(1,:)*[local_edges(2,2); -local_edges(2,1)];        
        %~ e_e = local_edges*local_edges';
        %~ C = [ eye(3,3), -e_e/2.0/S; zeros(3,3), [-1,-1,1;1,-1,-1;-1,1,-1]*diag(diag(e_e))/4.0/S ]; 

        %~ %the three functional value of each base function.
        %~ for domain_edge_ind=1:3
            %~ elt_edge_ind = is_node_on_edge( node(t,:), domain, domain_edge_ind); 
            %~ if elt_edge_ind>0
                 %~ if local_edges(elt_edge_ind,:)*domain_edges(domain_edge_ind,:)'>0
                    %~ two_edges_sign = edge_index_sign(elt_edge_ind+3);
                %~ else
                    %~ two_edges_sign = -1*edge_index_sign(elt_edge_ind+3);
                %~ end
                %~ g_index = [n_num+ abs(edge_index(elt_edge_ind) ) ];
                %~ OutM1( domain_edge_ind, g_index ) = two_edges_sign;               
            %~ end
        %~ end
 
        %~ local_tri=t;
        %~ clear t;
        %~ %the coefficients of three functions under area coordinates.
        
        %~ %coefficient for Global L1L2, L2L3,L3L1
        %~ st = st_coord(local_tri,:); s=st(:,1); t=st(:,2);
        %~ st_midpoint = 0.5*(st+st([2,3,1],:));
        %~ s_m=st_midpoint(:,1); t_m=st_midpoint(:,2);
        %~ v_G_L1L2_p = s.*t;
        %~ v_G_L2L3_p = (1-s-t).*t;
        %~ v_G_L3L1_p = s.*(1-s-t);
        %~ v_G_L1L2_m = s_m.*t_m  ;
        %~ v_G_L2L3_m = (1-s_m-t_m).*t_m;
        %~ v_G_L3L1_m = s_m.*(1-s_m-t_m);
        
        %~ coef_G_L1L2_m = 4*v_G_L1L2_m - 2*(v_G_L1L2_p + v_G_L1L2_p([2,3,1],:) );
        %~ coef_G_L2L3_m = 4*v_G_L2L3_m - 2*(v_G_L2L3_p + v_G_L2L3_p([2,3,1],:) );
        %~ coef_G_L3L1_m = 4*v_G_L3L1_m - 2*(v_G_L3L1_p + v_G_L3L1_p([2,3,1],:) );
        
        %~ %3 rows for 3 functions.
        %~ u = [v_G_L1L2_p,v_G_L2L3_p,v_G_L3L1_p;coef_G_L1L2_m,coef_G_L2L3_m,coef_G_L3L1_m]';
        
        %~ %the coefficients of three functions under local morley basisc.
        %~ u_morley = u*C*diag(edge_index_sign);
        
        %~ g_index = [local_tri, n_num+ abs(edge_index) ];
        
        %~ OutM2(:,g_index) = coeff_Function*u_morley;
 
    %~ end
    %~ H=OutM2';
%~ end
