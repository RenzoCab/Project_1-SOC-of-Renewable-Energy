function [nQ,nb,nc,nd,n_sets,n_bins] = Change_Index(Q,b,c,d,sets,bins)

    % We know that 3 is in sets and 4 is not. Then length(sets) > 0.

    nc = c;
    nQ = Q;
    nb = b;
    nd = d;
    nQ(1,1) = Q(3,3);
    nQ(1,2) = Q(3,4);
    nQ(2,1) = Q(4,3);
    nQ(3,3) = Q(1,1);
    nQ(3,4) = Q(1,2);
    nQ(4,3) = Q(2,1);
    nb(1) = b(3);
    nb(2) = b(4);
    nb(3) = b(1);
    nb(4) = b(2);
    nd(1) = d(3);
    nd(2) = d(4);
    nd(3) = d(1);
    nd(4) = d(2);
    
    if ismember(3,sets)
        n_sets(1) = 1;
        n_bins(1) = bins(sets==3);
    end
    if ismember(4,sets) % We know that 4 is not in sets.
        n_sets(end+1) = 2;
        n_bins(end+1) = bins(sets==4);
    end
    if ismember(1,sets)
        n_sets(end+1) = 3;
        n_bins(end+1) = bins(sets==1);
    end
    if ismember(2,sets)
        n_sets(end+1) = 4;
        n_bins(end+1) = bins(sets==2);
    end
    if ismember(5,sets)
        n_sets(end+1) = 5;
        n_bins(end+1) = bins(sets==5);
    end
    if ismember(6,sets)
        n_sets(end+1) = 6;
        n_bins(end+1) = bins(sets==6);
    end
    
end