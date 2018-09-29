function yes_or_no = is_edge_on_bd_edge_of_tri_domain(elt_edge, triangle, bd_edge_ind)
    
    %triangle: domain [x1,y1;x2,y2; x3,y3]
    %edge: format like [x1,y1;x2,y2].
    %bd_edge_ind: 1 or 2 or 3
    
    bd_edge = triangle( mod(bd_edge_ind,3)+1,:) - triangle( mod(bd_edge_ind-1,3)+1,:); 
    %~ if bd_edge_ind==2
        %~ '--'
        %~ elt_edge
        %~ bd_edge
    %~ end
    vecs=[elt_edge(:,1) - triangle( mod(bd_edge_ind,3)+1,1), elt_edge(:,2) - triangle( mod(bd_edge_ind,3)+1,2)];
    test = abs(vecs*[-bd_edge(1,2),bd_edge(1,1)]') < 1E-6;        
    if sum(test)>=2
        yes_or_no = 1;
        return; % the edge is on the boundary edge of triangle domain.
    end
    yes_or_no = 0;
end
