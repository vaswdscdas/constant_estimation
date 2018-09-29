function [OutA0, OutA1, OutA2,OutA_LM]=create_matrix_morley(tri,node,edge)
%% Inner product calcualtion of Fujino-Morley basis function (P2)

%  Sample value of parameters.    
%  tri = [1,2,3];
%  edge = [1,2;2,3;1,3];
%  node = [0,0;1,0;0,1];
%
%  Output matrix:
%  OutA0: \int_\Omega \phi_i \phi_j dxdy
%  OutA1: \int_\Omega \nabla \phi_i \nabla \phi_j dxdy
%  OutA2: \int_\Omega \nabla^2 \phi_i \nabla^2 \phi_j dxdy

%% Xuefeng LIU 
%% First version: 2012/09/21-2012/09/30

% Document for the matrix creating: <Dropbox>/Projects/FEM_Code_LIB/tex/main.pdf
% Detailed description for L2 inner product only.
% 2018/05/08
    
% Interval version: 
% If the INTERVAL_MODE variable is defined, then the output matrix will be interval ones.
% 2018/09/24

    e_num = size(edge,1); t_num = size(tri,1); n_num = size(node,1);
    tri_by_edge = zeros(t_num,3);

    edge = sort(edge,2);
    edge_value_index = edge*[n_num,1]';
    num = n_num + e_num;
    OutA0 = I_sparse( num, num ); OutA1 = I_sparse( num, num );  OutA2 = I_sparse( num, num );
    OutA_LM = I_sparse(num,num);
    
    for k=1:t_num  %create the table of element-edges relation.
        for l = 1:3  %the l-th edge.
            node_ind =[ mod(l,3)+1, mod(l+1,3)+1];
            local_value = sort( tri(k, node_ind) )*[n_num,1]';
            ind = find( edge_value_index - local_value == 0);
            tri_by_edge(k,l) = - ind * ( 2*( norm( tri(k,node_ind)-edge(ind,:) ) >0 ) - 1);
        end
    end
    for k=1:t_num
        t=tri(k,:);
        edge_index = tri_by_edge(k,:); %this can be negative.
        edge_index_sign = [1,1,1, (2*(edge_index>0) -1 )];
        edges = node( t([3,1,2]), : ) - node( t([2,3,1]), : );
        S =  edges(1,:)*[edges(2,2); -edges(2,1)]/2;
        e_e = edges*edges';

        %%  ------------   Method 1 to calculate P [START] ---------------------
        %%  This method needs inverse computation of a 6*6 matrix C.
        %C = [ eye(3,3), -e_e/2.0/S; zeros(3,3), [-1,-1,1;1,-1,-1;-1,1,-1]*diag(diag(e_e))/4.0/S ]; 
        %P = diag(edge_index_sign)*inv(C);
        %%  ------------   Method 1 to calculate P [END] ---------------------

        %%  ------------   Method 2 to calculate P [START] ---------------------
        TMP=diag( 1 ./ diag(e_e))*[1,0,1;1,1,0;0,1,1];
        CInv_22=-2*S*TMP;
        CInv_12=-e_e*TMP;
        CInv = [eye(3,3),CInv_12;zeros(3,3),CInv_22];
        %%  ------------   Method 2 to calculate P [END] ---------------------
        
        P = diag(edge_index_sign)*CInv;
                
        g_index = [t, n_num+ abs(edge_index) ];

        OutA_LM( g_index, g_index ) = OutA_LM( g_index, g_index ) + diag([S/3,S/3,S/3,0,0,0]);
        OutA0( g_index, g_index ) = OutA0( g_index, g_index ) + P*get_local_matrix_L2(S)*P';
        OutA1( g_index, g_index ) = OutA1( g_index, g_index ) + P*get_local_matrix_H1(e_e,S)*P';
        OutA2( g_index, g_index ) = OutA2( g_index, g_index ) + P*get_local_matrix_H2(e_e,S)*P';
    end
end

function local_matrix = get_local_matrix_L2(S)    %L_2 norm
    tmp = I_ones(3,3) + I_eye(3,3);
    local_matrix = [15*tmp,3*tmp; 3*tmp,tmp]*S/180;    %left-top
end

function local_matrix = get_local_matrix_H1(e_e,S) %H_1 norm
    local_matrix = I_zeros(6,6);
    local_matrix(1:3, 1:6) = [e_e, -e_e(:,[3,1,2])/3]/(4*S); %top 3 row
    local_matrix(4:6, 1:3) = local_matrix(1:3, 4:6)'; %left-bottom
    for k=1:3     %right-bottom    
        k1 = k; k2 = mod(k1,3)+1; k3 = mod(k2,3)+1;
        local_matrix( k1+3, k1+3) = ( e_e(k1,k1) - e_e(k3,k2) )/(24*S);
        local_matrix(k1+3,k2+3) = e_e(k1,k3)/(24*S);
        local_matrix(k2+3,k1+3) = local_matrix(k1+3,k2+3);
    end
end

function local_matrix = get_local_matrix_H2(e_e,S)  %H_2 norm
    local_matrix = I_zeros(6,6);
    for k=1:3
        k1 = k; k2 = mod(k1,3)+1; k3 = mod(k2,3)+1;
        local_matrix( k1+3, k1+3) = e_e(k1,k2)^2/(4*S^3) + 0.5/S;
        local_matrix( k1+3, k2+3) = e_e(k1,k2)*e_e(k2,k3)/(4*S^3) - 0.5/S;
        local_matrix( k2+3, k1+3) = local_matrix(k1+3, k2+3);
    end
end

% % % % 
% % % % function local_matrix = get_local_matrix_H1_complexform(edge,S)
% % % %     local_matrix=zeros(6,6);
% % % %     %left-top
% % % %     local_matrix(1:3, 1:3) = [ edge(1,:)*edge(1,:)', edge(1,:)*edge(2,:)', edge(1,:)*edge(3,:)';
% % % %                                edge(2,:)*edge(1,:)', edge(2,:)*edge(2,:)', edge(2,:)*edge(3,:)';
% % % %                                edge(3,:)*edge(1,:)', edge(3,:)*edge(2,:)', edge(3,:)*edge(3,:)';
% % % %                                 ]/4.0/S;
% % % % 
% % % %     %right-top
% % % %     local_matrix(1:3, 4:6) = -[ edge(1,:)*edge(3,:)', edge(1,:)*edge(1,:)', edge(1,:)*edge(2,:)';
% % % %                                edge(2,:)*edge(3,:)', edge(2,:)*edge(1,:)', edge(2,:)*edge(2,:)';
% % % %                                edge(3,:)*edge(3,:)', edge(3,:)*edge(1,:)', edge(3,:)*edge(2,:)';
% % % %                                 ]/12.0/S;
% % % %     %left-bottom
% % % %     local_matrix(4:6, 1:3) = local_matrix(1:3, 4:6)';
% % % %     
% % % %     %right-bottom    
% % % %     ind_list=[1,2,3;2,3,1;3,1,2];
% % % %     for k=1:3
% % % %         k1=ind_list(k,1);k2=ind_list(k,2);k3=ind_list(k,3);
% % % %         tmp = edge([k1,k2],:)*edge([k1,k2],:)'; 
% % % %         local_matrix( k1+3, k1+3) = (tmp(1,1)+tmp(2,2)+tmp(1,2))/24.0/S;
% % % %         local_matrix(k1+3,k2+3) = (edge(k1,:)*edge(k3,:)')/24.0/S;
% % % %         local_matrix(k2+3,k1+3) = local_matrix(k1+3,k2+3);
% % % %     end 
% % % %     
% % % % end
% % % % 
% % % % 
% % % % 
% % % % 
% % % % function local_matrix = get_local_matrix_H2_complateform(edge,S)
% % % %     local_matrix=zeros(6,6);
% % % %     r_edge=[ edge(:,2), -edge(:,1) ];
% % % %     e_dot_e = edge*edge'; 
% % % %     r_e_dot_e = r_edge*edge'; %roated e by 90 degree    
% % % %     ind_list=[1,2,3;2,3,1;3,1,2];
% % % %     for k=1:3
% % % %         k1=ind_list(k,1);k2=ind_list(k,2);k3=ind_list(k,3);
% % % %         for t=1:3
% % % %             t1=ind_list(t,1);t2=ind_list(t,2);t3=ind_list(t,3);
% % % %             local_matrix( k+3, t+3) = 2.0*edge(k1,2)*edge(k2,2)*edge(t1,2)*edge(t2,2);
% % % %             local_matrix( k+3, t+3) = local_matrix( k+3, t+3) + ( edge(k1,1)*edge(k2,2)+ edge(k1,2)*edge(k2,1))*(edge(t1,1)*edge(t2,2)+ edge(t1,2)*edge(t2,1)) ;
% % % %             local_matrix( k+3, t+3) = local_matrix( k+3, t+3) + 2.0*edge(k1,1)*edge(k2,1)*edge(t1,1)*edge(t2,1);
% % % %         end     
% % % %     end
% % % %     local_matrix = local_matrix /8.0/S^3;            
% % % % end
