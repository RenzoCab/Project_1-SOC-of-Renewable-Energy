function [newQ,newb,newc,newd] = Mat_Red(Q,b,c,d,i,x)

    % length(Q)=length(b)=length(d) is necessary.

    n = length(Q);
    
    for j=1:n
        if j>i
            newb(j-1) = b(j) + (Q(i,j)+Q(j,i))*x;
            newd(j-1) = d(j);
        elseif j ~= i
            newb(j) = b(j) + (Q(i,j)+Q(j,i))*x;
            newd(j) = d(j);
        end
        
        for k=1:n
            if j>i && k>i
                newQ(j-1,k-1) = Q(j,k);
            elseif j>i
                newQ(j-1,k) = Q(j,k);
            elseif k>i
                newQ(j,k-1) = Q(j,k);
            elseif j ~= i && k ~= i
                newQ(j,k) = Q(j,k);
            end
        end
    end

    newc = c + b(i)*x + Q(i,i)*x^2;

end