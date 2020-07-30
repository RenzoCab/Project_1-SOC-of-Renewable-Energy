function [nQ,nb,nc,nd] = Mat_Red_Bins(Q,b,c,d,sets,bins)

    if ~isempty(bins)

        for k = 1:length(bins)
            if k == 1
                [nQ,nb,nc,nd] = Mat_Red(Q,b,c,d,sets(1),bins(1));
            else
                [nQ,nb,nc,nd] = Mat_Red(nQ,nb,nc,nd,sets(k)+1-k,bins(k));
            end  
        end
    
    else
        
        nQ = Q;
        nb = b;
        nd = d;
        nc = c;
        
    end

end