function [nQ,nb,nc,nd] = App_Cond(Q,b,c,d,i,k,m)

    l = length(d);
    for j=1:l
        
        if j<m
            nd(j) = d(j);
        elseif j == m
            nd(j) = d(j)-d(j+1)*k(3)/k(2);
        elseif j>m+1
            nd(j-1) = d(j);
        end
        nd = nd';
        
        if j<m
            nb(j) = b(j);
        elseif j == m
            nb(j) = b(j)+Q(j+1,j)*2*(k(1)+i)/k(2);
        elseif j>m+1
            nb(j-1) = b(j);
        end
        nb = nb';
        
    end
    
    nQ = Q;
    nQ(:,m+1) = [];
    nQ(m+1,:) = [];
    nQ(m,m) = nQ(m,m)-Q(m+1,m)*2*k(3)/k(2);
    
    nc = c+d(5)*(i+k(1))/k(3);

end