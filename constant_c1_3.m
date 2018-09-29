%-----------------------------------------------------------------------
% Calcualte the Lagrange interpolation constant  over triangle T, which 
% has the middle points of edges as p12,p23,p31.
% 
% $\Pi u$ is a linear function satisfying
% (\PI_2 u) (p_i) = u(p_i), i=12,23,31
%
% The interpolation error estimation with constant C:
% || u - \PI u || <= C^0 |u|_2 for u in H^2(T)
% | u - \PI u |_1 <= C^1 |u|_2 for u in H^2(T)
%
% Notation : 
% C^0: C_1^{(3,0)}(T)
% C^1: C_1^{(3,1)}(T)
%-----------------------------------------------------------------------
%
% Requirement of mesh:
% The middile points of edges must be located on the grid of triangulation.
%
%
% Parameters of functions:
%
% mpath: path of mesh files.
% norm_idx: 0 or 1. (C^0： norm_idx=0, C^1: norm_idx=1）
%
% Xuefeng LIU 2012/10/05

% Check before publishing [in process]
% Started on 2018/05/08
% Interval version: 2018/09/29


function [c_value]=cal_constant_c1_3(mpath, norm_idx)

    tri = load([mpath,'tri_n.dat']); node = load([mpath,'node.dat']);
    edge = load([mpath,'edge.dat']);

    domain = load([mpath,'domain.dat']);
    domain_mid_p = 0.5*(domain + domain([2,3,1],:));

    ind1 = find( abs(node(:,1) - domain_mid_p(1,1) ) + abs(node(:,2) - domain_mid_p(1,2)) < 1E-10 );
    ind2 = find( abs(node(:,1) - domain_mid_p(2,1) ) + abs(node(:,2) - domain_mid_p(2,2)) < 1E-10 );
    ind3 = find( abs(node(:,1) - domain_mid_p(3,1) ) + abs(node(:,2) - domain_mid_p(3,2)) < 1E-10 );
    if length(ind1)*length(ind2)*length(ind3) < 1
        error('The middile points of edges must be located on the grid of triangulation.')
    end

    %The precision of node points, which will be used in rigourous computing.
    %This is only used for rigorous computing.    
    node_precision = 1E-14;
    node = I_set_interval_precision(node, node_precision);
    
    display(sprintf('node is interval?  %d \n', isintval(node)));
    
    h = get_max_edge_length(edge,node);    
    display(sprintf('mesh size: %s \n', I_sup(h)));

    node_num=size(node,1);
    [A0,A1,A2]=create_matrix_morley(tri,node,edge);

    n = size(A0,1);
    ind = 1:n;
    ind([ind1,ind2,ind3])=[];
    NA = A2(ind,ind);

    if norm_idx == 0
        NB = A0(ind,ind);
        Ch = (I_intval(0.0736)*h)^2;
    else
        NB = A1(ind,ind);
        Ch = (I_intval(0.1887)*h);

    end
  
    [v,lambda] = est_min_eig(NA,NB);
    lambda_lower = lambda/(1+(Ch)^2*lambda);   
    
    c_value =1./sqrt(lambda_lower);

end
