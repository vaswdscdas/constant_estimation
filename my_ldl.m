function [L,D, P] = my_ldl(A)
    if exist('ldl')
        [L,D,P]= ldl( mid(A) );
    else
        'error: TODO'
    end
end
