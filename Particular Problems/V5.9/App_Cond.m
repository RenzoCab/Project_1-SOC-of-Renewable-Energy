function [nQ,nb,nc,nd] = App_Cond(Q,b,c,d,i,k)

    l = length(d);
    for j=1:l
        
        if j<i
            nd(j) = d(j);
        elseif j == i
            nd(j) = d(j)-d(j+1)*k(3)/k(2);
        elseif j>i+1
            nd(j-1) = d(j);
        end
        nd = nd';
        
        if j<i
            nb(j) = b(j);
        elseif j == i
            nb(j) = b(j)+Q(j+1,j)*2*(k(1)+i)/k(2);
        elseif j>i+1
            nb(j-1) = b(j);
        end
        nb = nb';
        
    end
    
    nQ = Q;
    nQ(:,i+1) = [];
    nQ(i+1,:) = [];
    nQ(i,i) = nQ(i,i)-Q(i+1,i)*2*k(3)/k(2);
    
    nc = c;

end