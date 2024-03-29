close all;
clear all;
clc;

% Main Function for SOC.
% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email: renzo.caballerorosas@kaust.edu.sa caballerorenzo@hotmail.com
% Website: None.
% December 2019; Last revision: 30/12/2019.

Day = '20180302';
load([pwd '/Simulations/table_',Day,'.mat']);
lambda21 = saveTable(5,1:2);
lambda32 = saveTable(5,3:4);
expT = 6;
expS = 2;

% set(0,'DefaultFigureVisible','off');
% Optimal_Solution(Day,2,1,0,...
% 	1,1,1,expS,expT,1,1,...
%  	lambda21,lambda32,Day,1);

close all;
set(0,'DefaultFigureVisible','off');
Optimal_Solution(Day,4,1,0,...
	1,1,1,expS,expT,1,1,...
 	lambda21,lambda32,Day,1);

% close all;
% set(0,'DefaultFigureVisible','off');
% Optimal_Solution(Day,6,1,0,...
% 	1,1,1,expS,expT,1,1,...
%  	lambda21,lambda32,Day,1);

% close all;
% set(0,'DefaultFigureVisible','off');
% Optimal_Solution(Day,7,1,0,...
% 	1,1,1,expS,expT,1,1,...
%  	lambda21,lambda32,Day,1);