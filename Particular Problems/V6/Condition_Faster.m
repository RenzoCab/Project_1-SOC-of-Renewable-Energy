function [x4,x5,x8] = Condition_Faster(Sensibility,k,ToCheck,BinsToCheck)
    
    if ToCheck == 9
        
        % BinsToCheck = [x4,x5,0].
        x4 = BinsToCheck(1);
        x5 = BinsToCheck(2);
        x8 = x5*k(3) + x4*k(2) - k(1);
        
        if x8 > 1 || x8 < 0
            x4 = 2;
        end
        
    elseif ToCheck == 12
        
        % BinsToCheck = [x4,x8,0].
        x4 = BinsToCheck(1);
        x8 = BinsToCheck(2);
        x5 = (x8 - x4*k(2) + k(1)) / k(3);
        
        if x5 > 1 || x5 < 0
            x4 = 2;
        end
        
    elseif ToCheck == 13
        
        % BinsToCheck = [x5,x8,0].
        x5 = BinsToCheck(1);
        x8 = BinsToCheck(2);
        x4 = (x8 - x5*k(3) + k(1)) / k(2);
        
        if x4 > 1 || x4 < 0
            x4 = 2;
        end
        
    elseif ToCheck == 17
        
        % BinsToCheck = [x4,x5,x8].
        if abs(BinsToCheck(3)-BinsToCheck(2)*k(3)-BinsToCheck(1)*k(2)+k(1)) < Sensibility
            x4 = BinsToCheck(1);
            x5 = BinsToCheck(2);
            x8 = BinsToCheck(3);
        else
            x4 = 2;
            x5 = 0;
            x8 = 0;
        end

    end

end