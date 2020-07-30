function [newQ,newb,newc,newd] = Sys_Builder(Q,b,c,d,sets,bins)

    newQ = Q;
    newb = b;
    newc = c;
    newd = d;
    
    for i=1:length(bins)
        [newQ,newb,newc,newd] = Mat_Red(newQ,newb,newc,newd,sets(i)+1-i,bins(i));
    end
    
end