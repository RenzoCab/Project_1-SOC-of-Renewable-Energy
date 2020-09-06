function [A] = toRun(option)

%     set(0,'DefaultFigureVisible','off');
    set(0,'DefaultFigureVisible','on');
    
    %%%%%%%%%% 12/08/2019 %%%%%%%%%%

    dateOption1 = '20190103';
    discTimeOP1 = 6;
    discSpaceOP1 = 2;
    savePlotsOP1 = 0;
    
    if option == 1.3
        
        date = dateOption1;

        % This is to test the day 20180303. I will run this remotely.
        A = Optimal_Solution(['TEST_',date],3,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(A);
    
    elseif option == 1.2
        
        date = dateOption1;

        A = Optimal_Solution(['TEST_',date],2,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(A);
        
    elseif option == 1.25
        
        date = dateOption1;

        A = Optimal_Solution(['TEST_',date],1,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(A);
        B = Optimal_Solution(['TEST_',date],2,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(B);
        
        DG = abs(abs(A-B)/A);
        disp(DG);
        disp(sign(A*B));
        
    elseif option == 1.6
        
        date = dateOption1;
        discLambdaOption1 = 1;
        discLambdaOption1 = 999;
        
        % This is to test the day 20180303. I will run this remotely.
        A = Optimal_Solution(['TEST_',date],6,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,discLambdaOption1);
        disp(A);
    
    elseif option == 1.4
        
        date = dateOption1;

        % This is to test the day 20180303. I will run this remotely.
        Optimal_Solution(date,4,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,1);
    
    elseif option == 1.7
        
        date = dateOption1;

        % This is to test the day 20180303. I will run this remotely.
        Optimal_Solution(date,7,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,1);
    
    elseif option == 1.1
        
        date = dateOption1;

        % This is to test the day 20180303. I will run this remotely.
        A = Optimal_Solution(['TEST_',date],1,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(A);
    
    elseif option == 1.65
        
        date = dateOption1;

        A = Optimal_Solution(['TEST_',date],1,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(A);
        
        B = Optimal_Solution(['TEST_',date],6,1,0,...
        1,1,savePlotsOP1,discSpaceOP1,discTimeOP1,1,1,...
        0,0,date,999);
        disp(B);
        
        DG = abs(abs(A-B(1))/A);
        disp(DG);
        disp(sign(A*B(1)));
    
    end

    %%%%%%%%%% 13/08/2019 %%%%%%%%%%

    if option == 2

%         initilDate = '20180302';
        initilDate = '20200601';
        dateTime = datetime(initilDate,'InputFormat','yyyyMMdd');

        finalVector = [2020,10,1,0,0,0]; % Final date as a array.
%         finalDate = datetime(finalVector);
        finalDate = datestr(today-3,'yyyymmdd');
        finalDate = datetime(finalDate,'InputFormat','yyyyMMdd');
        errors = {};

        while dateTime <= finalDate
            
            timeString = datestr(dateTime,'yyyymmdd');
            if 0 == exist(['Historical/dailyData/Day_',timeString,'.mat'])
                
                try
                    date = datestr(dateTime,'yyyymmdd');
                    dateVec = datevec(dateTime);
                    dateVec = dateVec(1:3);
                    dataMaker(dateVec,dateVec);
%                     realPlotsMaker(dateVec,dateVec)
%                     Optimal_Solution(date,3,1,0,...
%                     1,1,1,2,8,1,1,...
%                     0,0,date,999);
                catch
                    errors{end+1} = dateTime;
                end
            
            else
                disp([timeString,' exists!']);
            end

        dateTime = dateTime + days(1);

        end

    end
    
    %%%%%%%%%% 20/08/2019 %%%%%%%%%%

    % The idea of this part is to simulate the two days ago (w.r.t. today)
    % data. We create the matrix and after simulate.
    
    if option == 3

        twoDaysAgo = datestr(today-3,'yyyymmdd');
        twoDaysAgoVec = datevec(today-3);
        twoDaysAgoVec = twoDaysAgoVec(1:3);
        errors = {};
        
        % I changed to 3 days ago.
        
        if not(isfile(['./Simulations/Simulation_',twoDaysAgo,'/101.eps']))
            try
                dataMaker(twoDaysAgoVec,twoDaysAgoVec);
                Optimal_Solution(twoDaysAgo,3,1,0,...
                1,1,1,2,8,1,1,...
                0,0,twoDaysAgo,999);
            catch
                disp('Error in toRun(3)!');
            end
        end
        
    end
    
    %%%%%%%%%% 21/09/2019 %%%%%%%%%%
    
    if option == 4
        
        discTimeOP4 = 8;
        discSpaceOP4 = 3;
        
        initilDate = '20190102';
        dateTime = datetime(initilDate,'InputFormat','yyyyMMdd');

        finalDate = '20200610';
        finalDate = datetime(finalDate,'InputFormat','yyyyMMdd');
        errors = {};

        while dateTime <= finalDate

            date = datestr(dateTime,'yyyymmdd');
            
            if not(isfile(['./Simulations/Simulation_',date,'/101.eps']))
                
                try
                    date = datestr(dateTime,'yyyymmdd');
                    
                    parpool('local',40);
                    
                    Optimal_Solution(date,1,1,0,...
                    1,1,1,discSpaceOP4,discTimeOP4,1,1,...
                    0,0,date,999);
                
                    delete(gcp('nocreate'));    
                
                catch
                    errors{end+1} = dateTime;
                    disp(['Error on ',date])
                end
                
            else
                
                disp([date,' is done!']);
            
            end

            dateTime = dateTime + days(1);

        end
        
    end
    
    %%%%%%%%%% 21/09/2019 %%%%%%%%%%

    dateOption5 = '20190102';
    discTimeOP1 = 5;
    discSpaceOP1 = 1;
    
    if option == 5

        for i = 0:2
%             Optimal_Solution(['Convergence_',num2str(i),'_',dateOption5],3,1,0,...
%                 1,1,1,discSpaceOP1+i,discTimeOP1+i,1,1,...
%                 0,0,dateOption5,999); % Simulation.
            A(i+1) = Optimal_Solution(['Convergence_',num2str(i),'_',dateOption5],2,1,0,...
                0,0,0,discSpaceOP1+i,discTimeOP1+i,1,1,...
                0,0,dateOption5,999); % Check output.
        end

    end
    
    %%%%%%%%%% 13/10/2019 %%%%%%%%%%
    
    if option == 6
        
        discTimeOP4 = 8;
        discSpaceOP4 = 3;
        
        initilDate = '20190102';
        dateTime = datetime(initilDate,'InputFormat','yyyyMMdd');

        finalDate = '20190228';
        finalDate = datetime(finalDate,'InputFormat','yyyyMMdd');
        errors = {};

        while dateTime <= finalDate

            try
                date = datestr(dateTime,'yyyymmdd');
                    Optimal_Solution(date,1,1,0,...
                    1,1,1,discSpaceOP4,discTimeOP4,1,1,...
                    0,0,date,999);
                date = datestr(dateTime,'yyyymmdd');
                    Optimal_Solution(date,1,1,0,...
                    1,1,1,discSpaceOP4,discTimeOP4,1,1,...
                    0,0,date,999);
            catch
                errors{end+1} = dateTime;
            end

        dateTime = dateTime + days(1);

        end
        
    end
    
    %%%%%%%%%% END %%%%%%%%%%
    
end