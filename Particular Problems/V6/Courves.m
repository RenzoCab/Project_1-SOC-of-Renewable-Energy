close all;
clear all;
clc;

%% R_040319

% Matrix = xlsread('R_040319.xlsx',1);
Matrix = xlsread('Real_January-February_2019.xlsx',1);

F = {'Fecha' 'Salto Grande' 'Bonete' 'Baygorria' 'Palmar' 'Eólica' 'Solar' 'Térmica' 'Biomasa' 'Imp.Arg' 'Imp.Br.Riv' ...
    'Imp.Br.Mel' 'Demanda' 'Exp_Intercon_BR_MELO' 'Exp_Intercon_BR_RIVERA' 'Exp_Intercon_ARG' 'Imp_Intercon_BR_MELO' 'Imp_Intercon_BR_RIVERA' 'Imp_Intercon_AG_Imp'};

Fecha = {'01/01/2019','02/01/2019','03/01/2019','04/01/2019','05/01/2019','06/01/2019','07/01/2019','08/01/2019','09/01/2019','10/01/2019',...
    '11/01/2019','12/01/2019','13/01/2019','14/01/2019','15/01/2019','16/01/2019','17/01/2019','18/01/2019','19/01/2019','20/01/2019',...
    '21/01/2019','22/01/2019','23/01/2019','24/01/2019','25/01/2019','26/01/2019','27/01/2019','28/01/2019','29/01/2019','30/01/2019',...
    '31/01/2019','01/02/2019','02/02/2019','03/02/2019','04/02/2019','05/02/2019','06/02/2019','07/02/2019','08/02/2019','09/02/2019','10/02/2019',...
    '11/02/2019','12/02/2019','13/02/2019','14/02/2019','15/02/2019','16/02/2019','17/02/2019','18/02/2019','19/02/2019','20/02/2019',...
    '21/02/2019','22/02/2019','23/02/2019','24/02/2019','25/02/2019','26/02/2019','27/02/2019','28/02/2019'};

T = (length(Matrix(:,1))-1)/6/24;
dt = 1/24/6;
time = 0:dt:1;

for i = 0:T-1
    figure('Position', [10 10 900 420]);
    hh = area(time,[Matrix(1+145*i-i:145+145*i-i,9),Matrix(1+145*i-i:145+145*i-i,7),Matrix(1+145*i-i:145+145*i-i,6),...
        Matrix(1+145*i-i:145+145*i-i,2),Matrix(1+145*i-i:145+145*i-i,3),Matrix(1+145*i-i:145+145*i-i,4),Matrix(1+145*i-i:145+145*i-i,5),Matrix(1+145*i-i:145+145*i-i,8)]);
    hold on;
    P = plot(time,Matrix(1+145*i-i:145+145*i-i,13),'k'); P.LineWidth = 1;
    P = plot(time,Matrix(1+145*i-i:145+145*i-i,13)+Matrix(1+145*i-i:145+145*i-i,14)+Matrix(1+145*i-i:145+145*i-i,15)+Matrix(1+145*i-i:145+145*i-i,16),'m');
    P.LineWidth = 1;
    grid minor;
    legend(F{9},F{7},F{6},F{2:5},F{8},F{13},'Demanda + Exp.','location','eastoutside');
    hh(1).FaceColor = [1 .5 0]; % Biomasa.
    hh(2).FaceColor = [1 1 .5]; % Solar.
    hh(3).FaceColor = [0 1 0]; % Wind.
    hh(4).FaceColor = [.5 .5 1]; % Bonete.
    hh(5).FaceColor = [.5 .8 1]; % Baygorria.
    hh(6).FaceColor = [.2 .2 1]; % Palmar.
    hh(7).FaceColor = [.8 .8 1]; % SG.
    hh(8).FaceColor = [1 0 0]; % Fuel.
    ylabel('MW');
    xlabel('Time');
    title(['Historical Production ',Fecha(i+1)]);
    mkdir('Historical');
%     saveas(gcf,[pwd '/Historical/',num2str(i+1)],'epsc');
end