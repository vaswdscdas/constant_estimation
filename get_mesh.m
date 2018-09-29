function [tris,nodes,mesh_edges,domain]=get_mesh(varargin)
    tri_node= varargin{1};
    h=varargin{2};
	 
    edges= tri_node([3,1,2],:) - tri_node([2,3,1],:);
    lens = sqrt( edges(:,1).^2 + edges(:,2).^2);
    [max_lens,ind]= max(lens);
    ind =  ind(1,1); %%index for node on the contrary side of maximal edge length.
    node_ind=[ ind, mod(ind,3)+1, mod(ind+1,3)+1];
    domain = tri_node( node_ind,:);
    edges= domain([3,1,2],:) - domain([2,3,1],:);

    tmp = (sum( lens(node_ind(2:3) ).^2) -lens(ind)^2 )/2.0/lens(node_ind(2))/lens(node_ind(3));
    max_angle = acos( tmp );  

       
    n = round(1/h);
    if mod(n,2) == 1 
       n = n+1;
    end
    h1= 1.0/n; %%along e3 direction
    h2= 1.0/n; %%along e2 direction
    vec1 = edges(3,:);
    vec2 = - edges(2,:);


    tris=[];
    mesh_edges=[];
    start = 0;
    ref_p = tri_node(ind,:);


%    if max_angle < 3.1415*2/3.0

    if max_angle < 3.1416/2.0
    
        nodes=[];
        for k=0:n
            x =h1*(0:n-k);
            shift_p = ( ref_p + k*h2*vec2 );
            tmp = x'*vec1 ;
            new_nodes =[ tmp(:,1)+shift_p(1), tmp(:,2)+shift_p(2) ];
            nodes=[nodes; new_nodes];
        end
    
        for k=1:n
            basic_tri1 = [1,2,n+2-(k-1)];
            basic_tri2 = [2,n+2-(k-1)+1,n+2-(k-1)];
           for l=1:n - (k-1)
               mesh_edges = [mesh_edges; basic_tri1(1:2)+ l-1+start ]; 
               mesh_edges = [mesh_edges; basic_tri1(2:3)+ l-1+start ]; 
               mesh_edges = [mesh_edges; basic_tri1([1,3])+ l-1+start ]; 

                tris = [tris; basic_tri1+l-1+start];
            end
            for l=1:n-1- (k-1)
                tris = [tris; basic_tri2+l-1+start];
            end
            start = start + (n+1) - (k-1);
        end
    else
    
        nodes=[];
        for k=0:n
            x =h1*(0:n-k);
            shift_p = ( ref_p + k*h2*vec2 );
            tmp = x'*vec1;
            new_nodes =[ tmp(:,1)+shift_p(1), tmp(:,2)+shift_p(2) ];
            nodes=[nodes; new_nodes];
            if k<n 
                new_point_on_max_edge = (n-k-0.5)*h1*vec1 + h2*vec2*(k+0.5) + ref_p;
                nodes=[nodes; new_point_on_max_edge];
            end
        end

    
        for k=1:n
        
            basic_tri1 = [1,n+2-(k-1)+2, n+2-(k-1)+1 ];
            basic_tri2 = [1,2,n+2-(k-1)+2];
           for l=1:n -1 - (k-1)
                tris = [tris; basic_tri1+l-1+start];
                tris = [tris; basic_tri2+l-1+start];
                mesh_edges = [mesh_edges; [1,2] + l-1+start ]; 
                mesh_edges = [mesh_edges; [1,n+2-(k-1)+2] + l-1+start ]; 
                mesh_edges = [mesh_edges; [1,n+2-(k-1)+1] + l-1+start ]; 
            end
            tmp1 = [n-(k-1),n+2-(k-1), 2*n-2*k+4]+start;
            tmp2 = [n-(k-1),n-(k-1)+1,n+2-(k-1)]+start;            
            mesh_edges = [mesh_edges; tmp2([1,2]) ]; 
            mesh_edges = [mesh_edges; tmp2([2,3]) ]; 
            mesh_edges = [mesh_edges; tmp2([1,3]) ]; 
            mesh_edges = [mesh_edges; tmp1([2,3]) ]; 
            mesh_edges = [mesh_edges; tmp1([1,3]) ]; 
            
            tris = [tris; tmp1;tmp2 ];
            start = start + (n+1) - (k-1)+1;
        end
    end
 
    if length(varargin)==3
	    folder=varargin{3};
	    if folder(length(folder)) ~= filesep() 
		folder = folder+filesep();
	    end
	    save([folder,'tri_n.dat'],'tris','-ascii');
	    save([folder,'node.dat'],'nodes','-ascii', '-double');
	    save([folder,'edge.dat'],'mesh_edges','-ascii');
	    save([folder,'domain.dat'],'domain','-ascii');
    end
    

%    vertex =nodes;
%    face=[tris, tris(:,1)];
%	xvf = matrix(vertex(face,1),size(face,1),length(vertex(face,1))/size(face,1))';
%	yvf = matrix(vertex(face,2),size(face,1),length(vertex(face,1))/size(face,1))';
%    plot(xvf, yvf,'b-');
%    
    end

