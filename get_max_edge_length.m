function h = get_max_edge_length(edge,node)
% Given edge list and node list, find the largest ede length.
%
% Histrory
% Start: 2018/05/08, Xuefeng LIU 
% Bug fixed for vec calculation. 2018/09

 vec = node(edge(:,1),:) - node(edge(:,2),:);
 h = sqrt( max(sum( vec.^2, 2)) );
 
end
