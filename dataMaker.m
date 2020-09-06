function [] = dataMaker(initialVector,finalVector)

%     initialVector = [2018,3,2]; % Initial date as a array.
%     initialVector = [2019,1,1]; % Initial date as a array.
    initialVector = [initialVector,0,0,0];
    initialDate = datetime(initialVector); % Initial date as a Time variable.

%     finalVector = [2019,7,1]; % Final date as a array.
    finalVector = [finalVector,0,0,0];
    finalDate = datetime(finalVector); % Final date as a Time variable.

    date = initialDate;

    while date <= finalDate

        try
            dateVector = datevec(date);
            New_Real_Data(dateVector);
            disp(['Day ',datestr(date,'dd/mm/yyyy'),' was successful.'])
        catch
            warning(['There is something wrong with the date ',datestr(date,'dd/mm/yyyy')]);
        end

        date = date + days(1);
    
    end
    
end