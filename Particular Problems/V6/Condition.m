function [x4,x5,x8] = Condition(sets,bins,k)

    Sensibility = 0.01;

    if ismember(4,sets)&&(ismember(5,sets))&&(ismember(8,sets))
        
        if abs(bins(sets==8)-bins(sets==5)*k(3)-bins(sets==4)*k(2)+k(1)) < Sensibility
            x4 = bins(sets==4);
            x5 = bins(sets==5);
            x8 = bins(sets==8);
        else
            x4 = 2; % We use 2 because it is greater than 1.
            x5 = 2;
            x8 = 2;
        end
        
    elseif ismember(4,sets)&&(ismember(5,sets))
        x4 = bins(sets==4);
        x5 = bins(sets==5);
        x8 = x5*k(3) + x4*k(2) - k(1);
    elseif ismember(4,sets)&&(ismember(8,sets))
        x4 = bins(sets==4);
        x8 = bins(sets==8);
        x5 = (x8 - x4*k(2) + k(1)) / k(3);
    else % elseif ismember(5,sets)&&(ismember(8,sets))
        x5 = bins(sets==5);
        x8 = bins(sets==8);
        x4 = (x8 - x5*k(3) + k(1)) / k(2);
    end

    if sum([x4,x5,x8]>1) + sum([x4,x5,x8]<0)
        x4 = 1i;
        x5 = 1i;
        x8 = 1i;
    end

end