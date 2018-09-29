function [elt_edge_ind]=is_node_on_edge(element, triangle, edge_ind)
    %Find the edge index of element which is located on the specified edge (by edge_ind) of triangle.
    
    edge = triangle( mod(edge_ind+1,3)+1,:) - triangle( mod(edge_ind,3)+1,:) ; 
    
    vecs=[element(:,1) - triangle( mod(edge_ind,3)+1,1), element(:,2) - triangle( mod(edge_ind,3)+1,2)];
    
    test = abs(vecs*[-edge(1,2),edge(1,1)]')>1E-6;        
    if sum(test)>=2 
        elt_edge_ind = -1;
        return; % no edge of element on triangle.
    end
    [ind]=find( test >0 );
    elt_edge_ind = ind;
end
