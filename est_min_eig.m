% Given matrix A and B, find the minimum eigenalue of 
% Ax = lambda Bx
%
%
% Histrory
% Start: 2018/05/08, Xuefeng LIU 
% 

function v = est_min_eig(A,B)

   global INTERVAL_MODE;

   if INTERVAL_MODE
     B = hull(B, B');
     A = hull(A, A');

     [v,idx] = veigs( A, B, 'sm',1);
   else
     v = eigs( A, B, 1, 'sm');
   end

end
