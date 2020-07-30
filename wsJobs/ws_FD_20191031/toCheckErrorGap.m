close all;
clear all;
clc;

set(0,'DefaultFigureVisible','on');
doWePlot = 1;

Day = '20180302';
numDays = 600;
aux = {};
vectorGAP = zeros(1,numDays);
finalCostHJB = zeros(1,numDays);
finalCostToGo = zeros(1,numDays);
finalCostToGoReal = zeros(1,numDays);
kh = zeros(4,numDays);
In = zeros(4,numDays);
badDays{600} = {};

if 0 == exist([pwd '/PostPostProcessing/Plots'],'dir')
	mkdir([pwd '/PostPostProcessing/Plots']);
end

for i = 1:numDays
    try
        load([pwd '/','Simulations/dataFrom_6_',Day,'.mat']);
        load([pwd '/','Simulations/dataFrom_7_',Day,'.mat']);
        finalCostHJB(i) = dataFrom_6{1};
        finalCostToGo(i) = dataFrom_6{2};
        finalCostToGoReal(i) = dataFrom_7{1};
        kh(1,i) = dataFrom_6{3}; kh(2,i) = dataFrom_6{4};
        kh(3,i) = dataFrom_6{5}; kh(4,i) = dataFrom_6{6};
        In(1,i) = dataFrom_6{7}; In(2,i) = dataFrom_6{8};
        In(3,i) = dataFrom_6{9}; In(4,i) = dataFrom_6{10};
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

for i = numDays:-1:1 % Remove uncomputed days.
    if vectorGAP(i) == 0
        vectorGAP(i) = [];
        kh(:,i) = [];
        In(:,i) = [];
        badDays(i) = [];
        finalCostHJB(i) = [];
        finalCostToGo(i) = [];
        finalCostToGoReal(i) = [];
    end
end
xAxesVectorGap = 1:1:length(vectorGAP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Above we loaded the data.

figure(1);
Y = area(xAxesVectorGap,[ones(1,length(xAxesVectorGap))*3;...
    ones(1,length(xAxesVectorGap))*2;ones(1,length(xAxesVectorGap))*95]');
Y(1).FaceColor = [160 255 160]/255;
Y(2).FaceColor = [255 255 160]/255;
Y(3).FaceColor = [255 160 160]/255;
hold on;
plot(xAxesVectorGap,vectorGAP,'o-');
xlabel('Days');
ylabel('Relative Error (%)');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
ylim([0,100]);
box;

badDaysCounter = 1;
for i = 1:length(vectorGAP)
    if vectorGAP(i) > 5
%         text(xAxesVectorGap(i),vectorGAP(i),badDays{i});
        setBadDays{badDaysCounter} = badDays{i};
        badDaysCounter = badDaysCounter + 1;
    end
end

errorPorc = length(setBadDays)/length(vectorGAP)*100;
disp([num2str(errorPorc),'% of the days have a big error GAP.']);
title(['Discretization Relative Error (',num2str(errorPorc),'% are red)']);

figure(2);
histogram(finalCostHJB,200);
xlabel('Hundred Thousand USD (HJB)');
title('HJB initial value histogram');
grid minor;
box;

figure(8);
histogram(finalCostToGoReal,200);
xlabel('Hundred Thousand USD (HJB)');
title('Real Running Cost');
grid minor;
box;

figure(3);
histogram(abs(finalCostHJB),200);
xlabel('Hundred Thousand USD (HJB)');
title('Absolute HJB initial value histogram');
grid minor;
box;

figure(4);
hold on;
plot(xAxesVectorGap,kh(1,:),'LineWidth',3);
plot(xAxesVectorGap,kh(2,:),'LineWidth',1.5);
P = plot(xAxesVectorGap,kh(3,:),'--');
P.LineWidth = 1.5;
plot(xAxesVectorGap,kh(4,:),'LineWidth',1.5);
grid minor;
title('Water Value');
xlabel('Days');
ylabel('USD/MWh');
legend('Bonete','Baygorria','Palmar','SG');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
box;

figure(5);
hold on;
plot(xAxesVectorGap,In(1,:));
plot(xAxesVectorGap,In(2,:));
plot(xAxesVectorGap,In(3,:));
plot(xAxesVectorGap,In(4,:));
grid minor;
title('Natural Inflow');
xlabel('Days');
ylabel('$m^3/s$','Interpreter','latex');
legend('Bonete','Baygorria','Palmar','SG');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
box;

figure(6);
hold on;
plot(xAxesVectorGap,finalCostHJB);
plot(xAxesVectorGap,finalCostToGo);
plot(xAxesVectorGap,zeros(1,length(xAxesVectorGap)),'k--');
grid minor;
title('Value Function and Running Cost');
xlabel('Days');
ylabel('Hundred Thousand USD');
legend('Value Function','Running Cost');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
box;

figure(9);
hold on;
plot(xAxesVectorGap,finalCostToGo);
plot(xAxesVectorGap,finalCostToGoReal);
plot(xAxesVectorGap,zeros(1,length(xAxesVectorGap)),'k--');
grid minor;
title('Simulated and Historical Running Costs');
xlabel('Days');
ylabel('Hundred Thousand USD');
legend('Running Cost','Historical Running Cost');
xlim([xAxesVectorGap(1),xAxesVectorGap(end)]);
box;

auxFinalCostHJB = abs(finalCostHJB); % We remove the values greater than 5.
for i = length(auxFinalCostHJB):-1:1
    if auxFinalCostHJB(i) >= 5
        auxFinalCostHJB(i) = [];
    end
end
figure(7);
histogram(auxFinalCostHJB,200);
xlabel('Hundred Thousand USD (HJB)');
title('Absolute HJB (only [values < 5])');
grid minor;
box;

j = 1;
for i = length(vectorGAP):-1:1
    if vectorGAP(i) >= 5
        if finalCostHJB(i) * finalCostToGo(i) > 0 % I verify the same sign.
            absoluteError(j) = abs(finalCostHJB(i)-finalCostToGo(i));
            relativeError(j) = 100 * abs(finalCostHJB(i)-finalCostToGo(i)) / abs(finalCostHJB(i));
            j = j + 1;
        end
    end
end
figure(10);
plot(1:1:length(absoluteError),absoluteError,'o-');
hold on;
plot(1:1:length(absoluteError),0.2*ones(1,length(absoluteError)),'k--');
grid minor;
title('Absolute Error');
xlabel('Days');
ylabel('Hundred Thousand USD');
xlim([1,length(absoluteError)]);
box on;
figure(11);
plot(1:1:length(relativeError),relativeError,'o-');
grid minor;
title('Relative Error');
xlabel('Days');
ylabel('Relative Error (%)')
xlim([1,length(relativeError)]);
box on;

if doWePlot == 1
    for i = 1:11
        set(0,'CurrentFigure',i);
        saveas(gcf,[pwd '/PostPostProcessing/Plots/Plot_',num2str(i)],'epsc');
    end
end

save('setBadDays.mat','setBadDays');