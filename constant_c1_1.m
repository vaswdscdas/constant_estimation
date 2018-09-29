%-----------------------------------------------------------------------
% Calcualte the Lagrange interpolation constant over triangle T, which 
% has the vertices as p1,p2,p3.
% 
% $\Pi u$ is a linear function satisfying
% (\PI_2 u) (p_i) = u(p_i), i=1,2,3
%
% The interpolation error estimation with constant C:
% || u - \PI u || <= C^0 |u|_2 for u in H^2(T)
% | u - \PI u |_1 <= C^1 |u|_2 for u in H^2(T)
%
% Notation : 
% C^0: C_1^{(1,0)}(T)
% C^1: C_1^{(1,1)}(T)
%
%-----------------------------------------------------------------------
%
% Parameters of functions:
%
% mpath: path of mesh files.
% norm_idx: 0 or 1. (C^0： norm_idx=0, C^1: norm_idx=1）

%
% Xuefeng LIU 
% First version: 2014/09/9
% Check before publishing [in process]
% Started on 2018/05/08

function [c_value]=constant_c1_1(mpath, norm_idx)

    tri = load([mpath,'tri_n.dat']); node = load([mpath,'node.dat']);
    edge = load([mpath,'edge.dat']);

    %The precision of node points, which will be used in rigourous computing.
    %This is only used for rigorous computing.    
    node_precision = 1E-14;
    node = I_set_interval_precision(node, node_precision);
    
    display(sprintf('Is node in interval mode?  %d \n', isintval(node)));
    
    h = get_max_edge_length(edge,node);    
    display(sprintf('Mesh size: %s \n', I_sup(h)));
    
    domain = load([mpath,'domain.dat']);

    ind1 = find( ( abs(node(:,1)-domain(1,1) ) + abs(node(:,2)-domain(1,2)) ) < 1E-10 );
    ind2 = find( ( abs(node(:,1)-domain(2,1) ) + abs(node(:,2)-domain(2,2)) ) < 1E-10 );
    ind3 = find( ( abs(node(:,1)-domain(3,1) ) + abs(node(:,2)-domain(3,2)) ) < 1E-10 );
    node_num = size(node,1);

    [A0,A1,A2] = create_matrix_morley(tri,node,edge);
        
    n = size(A0,1); % DOF of Morley FEM space.
    ind = 1:n;
    ind([ind1,ind2,ind3]) = [];
    NA = A2(ind,ind);

    if norm_idx == 0
        NB = A0(ind,ind);
        Ch = (I_intval(0.0736)*h)^2;
    else
        NB = A1(ind,ind);
        Ch = (I_intval(0.1887)*h);
    end
   
    lambda = est_min_eig(NA,NB);
    lambda_lower = lambda/(1 + Ch^2*lambda);
    
    c_value = 1./sqrt(lambda_lower);

    return 
    
    % Eigenvector drawing.
    
    figure(2)
    hold off
    vec=zeros(node_num,1); vec(ind)=v(:,1);
    trisurf(tri,node(:,1), node(:,2),vec(1:node_num));

end

