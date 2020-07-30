function [values] = toRun(input)

% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email: renzo.caballerorosas@kaust.edu.sa caballerorenzo@hotmail.com
% Website: None.
% September 2019; Last revision: 23/09/2019.    
    
    if input == 1
        initilDate = '20190101';
        finalDate = '20190131';
    elseif input == 2
        initilDate = '20190201';
        finalDate = '20190228';
    elseif input == 3
        initilDate = '20190101';
        finalDate = '20190228';
    end
    
    dateTime = datetime(initilDate,'InputFormat','yyyyMMdd');
    finalDate = datetime(finalDate,'InputFormat','yyyyMMdd');
    values = {};

    while dateTime <= finalDate

        try
            date = datestr(dateTime,'yyyymmdd');
            load(['output_',date,'.mat']);
            for i = 1:2
                lambda21(i) = output{1}(i);
                lambda32(i) = output{1}(i+2);
            end
            
%             A = Optimal_Solution([date,'_ZeroLambdas'],7,1,0,...
%                 1,1,1,2,8,1,1,...
%                 0,0,date,999);
%             B = Optimal_Solution([date,'_OptimalLambdas'],2,1,0,...
%                 1,1,1,2,8,1,1,...
%                 lambda21,lambda32,date,999);
%             C = Optimal_Solution([date,'_ZeroLambdas'],3,1,0,...
%                 1,1,1,2,8,1,1,...
%                 0,0,date,999);
%             values{end+1} = [A,B,C];
            
%             A = Optimal_Solution([date,'_ZeroLambdas'],7,1,0,...
%                 1,1,1,2,8,1,1,...
%                 0,0,date,999); % Remove this and uncomment all.
            B = Optimal_Solution([date,'_OptimalLambdas'],2,1,0,...
                1,1,1,2,8,1,1,...
                lambda21,lambda32,date,999); % Remove this and uncomment all.
            
        catch
            values{end+1} = 'error';
        end

    dateTime = dateTime + days(1);

    end
    
%     save('costs.mat','values');
    
end