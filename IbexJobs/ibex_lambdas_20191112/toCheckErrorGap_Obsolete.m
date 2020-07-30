close all;
clear all;
clc;

set(0,'DefaultFigureVisible','off');

Day = '20180302';
numDays = 600;
aux = {};
vectorGAP = zeros(1,numDays);
xAxesVectorGap = 1:1:numDays;
finalCostHJB = zeros(1,numDays);
finalCostToGo = zeros(1,numDays);
badDays{600} = {};

for i = 1:numDays
    try
        
        Optimal_Solution(Day,6,1,0,...
        1,1,0,3,8,1,1,...
        0,0,Day,999);
    
        set(0,'CurrentFigure',102);
        fig = gcf;
        axObjs = fig.Children;
        axObjs.Children;
        aux{1} = ans;
        finalCostHJB(i) = aux{1}(3).YData(end);
        finalCostToGo(i) = aux{1}(2).YData(end);
        errorGAP = 100 * abs(finalCostHJB(i)-finalCostToGo(i))/abs(finalCostHJB(i));
        vectorGAP(i) = errorGAP;
        if errorGAP > 5
            badDays{i} = Day;
        end
    
    end
    matFormat = datetime(Day,'InputFormat','yyyyMMdd');
    matFormat = matFormat + days(1);
    Day = datestr(matFormat,'yyyymmdd');
end

continuacion = 1;
valBackward = numDays;
while continuacion
    if vectorGAP(valBackward) == 0
        vectorGAP(valBackward) = [];
        xAxesVectorGap(valBackward) = [];
        badDays(valBackward) = [];
        finalCostHJB(valBackward) = [];
        finalCostToGo(valBackward) = [];
        valBackward = valBackward - 1;
    else
        continuacion = 0;
    end
end

infoToPlot = {xAxesVectorGap,vectorGAP,badDays,finalCostHJB,finalCostToGo};
save('infoToPlot.mat','infoToPlot');

figure;
Y = area(xAxesVectorGap,[ones(1,length(xAxesVectorGap))*3;...
    ones(1,length(xAxesVectorGap))*2;ones(1,length(xAxesVectorGap))*95]');
Y(1).FaceColor = [160 255 160]/255;
Y(2).FaceColor = [255 255 160]/255;
Y(3).FaceColor = [255 160 160]/255;
hold on;
plot(xAxesVectorGap,vectorGAP,'o-');
xlabel('Days');
ylabel('Relative Error');
title('Discretization Relative Error');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
ylim([0,max(vectorGAP)]);

for i = 1:length(vectorGAP)
    if vectorGAP(i) > 5
        text(xAxesVectorGap(i),vectorGAP(i),badDays{i});
    end
end

figure;
plot(xAxesVectorGap,vectorGAP,'o-');
xlabel('Days');
ylabel('Relative Error');
title('Discretization Relative Error');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
ylim([0,max(vectorGAP)]);