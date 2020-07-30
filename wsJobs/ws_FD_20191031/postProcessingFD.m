close all;
clear all;
clc;

Day = '20180302';
count = 0;
badCount = 0;
badCount_2 = 0;
badCount_3 = 0;
badCount_R = 0;
badCount_2_R = 0;
badCount_3_R = 0;
badCount_L = 0;
badCount_2_L = 0;
badCount_3_L = 0;
badCount_A = 0;
badCount_2_A = 0;
badCount_3_A = 0;

for i = 1:600
    
    if 0 ~= exist([pwd '/Simulations/FD_',Day,'.mat'])
        load([pwd '/Simulations/FD_',Day,'.mat']);
        for j = 1:10
            auxVec(j) = allData{j}(1);
        end
        count = count + 4;
        auxVec(10) = auxVec(10)*2; % For central FD and normalization.
        disp(['Day: ',Day]);
        disp(['Lambda21_FD (1): ',num2str((auxVec(8)-auxVec(9))/auxVec(10)),...
            ' Theorem: ',num2str(allData{1}(2))]);
        disp(['Lambda21_FD (2): ',num2str((auxVec(6)-auxVec(7))/auxVec(10)),...
            ' Theorem: ',num2str(allData{1}(3))]);
        disp(['Lambda32_FD (1): ',num2str((auxVec(4)-auxVec(5))/auxVec(10)),...
            ' Theorem: ',num2str(allData{1}(4))]);
        disp(['Lambda32_FD (2): ',num2str((auxVec(2)-auxVec(3))/auxVec(10)),...
            ' Theorem: ',num2str(allData{1}(5))]);
        
        r1 = abs(((auxVec(8)-auxVec(9))/auxVec(10)-allData{1}(2))/allData{1}(2));
        r2 = abs(((auxVec(6)-auxVec(7))/auxVec(10)-allData{1}(3))/allData{1}(3));
        r3 = abs(((auxVec(4)-auxVec(5))/auxVec(10)-allData{1}(4))/allData{1}(4));
        r4 = abs(((auxVec(2)-auxVec(3))/auxVec(10)-allData{1}(5))/allData{1}(5));
        
        auxVec(10) = auxVec(10)/2;
        
        r1_R = abs(((auxVec(8)-auxVec(1))/auxVec(10)-allData{1}(2))/allData{1}(2));
        r2_R = abs(((auxVec(6)-auxVec(1))/auxVec(10)-allData{1}(3))/allData{1}(3));
        r3_R = abs(((auxVec(4)-auxVec(1))/auxVec(10)-allData{1}(4))/allData{1}(4));
        r4_R = abs(((auxVec(2)-auxVec(1))/auxVec(10)-allData{1}(5))/allData{1}(5));
        
        r1_L = abs(((auxVec(1)-auxVec(9))/auxVec(10)-allData{1}(2))/allData{1}(2));
        r2_L = abs(((auxVec(1)-auxVec(7))/auxVec(10)-allData{1}(3))/allData{1}(3));
        r3_L = abs(((auxVec(1)-auxVec(5))/auxVec(10)-allData{1}(4))/allData{1}(4));
        r4_L = abs(((auxVec(1)-auxVec(3))/auxVec(10)-allData{1}(5))/allData{1}(5));
        
        r1_A = min(r1_R,r1_L);
        r2_A = min(r2_R,r2_L);
        r3_A = min(r3_R,r3_L);
        r4_A = min(r4_R,r4_L);
        
        threshold = 0.03;
        threshold_2 = 0.05;
        threshold_3 = 0.1;
        badCount = badCount + (r1 >= threshold) + (r2 >= threshold) + (r3 >= threshold) + (r4 >= threshold);
        badCount_2 = badCount_2 + (r1 >= threshold_2) + (r2 >= threshold_2) + (r3 >= threshold_2) + (r4 >= threshold_2);
        badCount_3 = badCount_3 + (r1 >= threshold_3) + (r2 >= threshold_3) + (r3 >= threshold_3) + (r4 >= threshold_3);
        disp(['Relative errors: ',num2str(r1),', ',num2str(r2),', ',num2str(r3),', ',num2str(r4)]);
        
        badCount_R   = badCount_R   + (r1_R >= threshold)   + (r2_R >= threshold)   + (r3_R >= threshold)   + (r4_R >= threshold);
        badCount_2_R = badCount_2_R + (r1_R >= threshold_2) + (r2_R >= threshold_2) + (r3_R >= threshold_2) + (r4_R >= threshold_2);
        badCount_3_R = badCount_3_R + (r1_R >= threshold_3) + (r2_R >= threshold_3) + (r3_R >= threshold_3) + (r4_R >= threshold_3);
        
        badCount_L   = badCount_L   + (r1_L >= threshold)   + (r2_L >= threshold)   + (r3_L >= threshold)   + (r4_L >= threshold);
        badCount_2_L = badCount_2_L + (r1_L >= threshold_2) + (r2_L >= threshold_2) + (r3_L >= threshold_2) + (r4_L >= threshold_2);
        badCount_3_L = badCount_3_L + (r1_L >= threshold_3) + (r2_L >= threshold_3) + (r3_L >= threshold_3) + (r4_L >= threshold_3);
        
        badCount_A   = badCount_A   + (r1_A >= threshold)   + (r2_A >= threshold)   + (r3_A >= threshold)   + (r4_A >= threshold);
        badCount_2_A = badCount_2_A + (r1_A >= threshold_2) + (r2_A >= threshold_2) + (r3_A >= threshold_2) + (r4_A >= threshold_2);
        badCount_3_A = badCount_3_A + (r1_A >= threshold_3) + (r2_A >= threshold_3) + (r3_A >= threshold_3) + (r4_A >= threshold_3);
    end
    
    matFormat = datetime(Day,'InputFormat','yyyyMMdd');
    matFormat = matFormat + days(1);
    Day = datestr(matFormat,'yyyymmdd');
    
end
disp('---------------------------------------')
disp(['We processed ',num2str(count),' subgradients']);
disp('---------------------------------------')
disp(['Central FD:']);
disp([num2str((badCount/count)*100),'% have relative error grather than ',num2str(threshold*100),'%']);
disp([num2str((badCount_2/count)*100),'% have relative error grather than ',num2str(threshold_2*100),'%']);
disp([num2str((badCount_3/count)*100),'% have relative error grather than ',num2str(threshold_3*100),'%']);
disp('---------------------------------------')
disp(['Right-hand FD:']);
disp([num2str((badCount_R/count)*100),'% have relative error grather than ',num2str(threshold*100),'%']);
disp([num2str((badCount_2_R/count)*100),'% have relative error grather than ',num2str(threshold_2*100),'%']);
disp([num2str((badCount_3_R/count)*100),'% have relative error grather than ',num2str(threshold_3*100),'%']);
disp('---------------------------------------')
disp(['Left-hand FD:']);
disp([num2str((badCount_L/count)*100),'% have relative error grather than ',num2str(threshold*100),'%']);
disp([num2str((badCount_2_L/count)*100),'% have relative error grather than ',num2str(threshold_2*100),'%']);
disp([num2str((badCount_3_L/count)*100),'% have relative error grather than ',num2str(threshold_3*100),'%']);
disp('---------------------------------------')
disp(['Min(L,R) FD:']);
disp([num2str((badCount_A/count)*100),'% have relative error grather than ',num2str(threshold*100),'%']);
disp([num2str((badCount_2_A/count)*100),'% have relative error grather than ',num2str(threshold_2*100),'%']);
disp([num2str((badCount_3_A/count)*100),'% have relative error grather than ',num2str(threshold_3*100),'%']);