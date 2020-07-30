function [] = realPlotsMaker(initialVector,finalVector)

%     initialVector = [2016,11,3]; % Initial date as a array.
    initialVector = [initialVector,0,0,0];
    initialDate = datetime(initialVector); % Initial date as a Time variable.

%     finalVector = [2019,7,1]; % Final date as a array.
    finalVector = [finalVector,0,0,0];
    finalDate = datetime(finalVector); % Final date as a Time variable.

    date = initialDate;

    while date <= finalDate

        try
            dateString = datestr(date,'yyyymmdd');
            realBalancePlot(dateString); % The input must be a string yyyyMMdd.
            realBalancePlotCont(dateString); % The input must be a string yyyyMMdd.
            print(['Day ',datestr(date,'dd/MM/yyyy'),' was successful.'])
        catch
            warning(['There is something wrong with the date ',datestr(date,'dd/mm/yyyy')]);
        end

        date = date + days(1);
    
    end
    
end