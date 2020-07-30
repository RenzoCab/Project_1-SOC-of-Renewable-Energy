% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 14/11/2019

close all;
clear all;
clc;

extremes = 10;
dx = 0.01;

[X,Y] = meshgrid(-extremes:dx:extremes,-extremes:dx:extremes);
for i = 1:2*extremes/dx+1
    for j = 1:2*extremes/dx+1
        [f,g] = weirdCone([X(i,j),Y(i,j)]);
        Z(i,j) = f;
        G1(i,j) = g(1);
        G2(i,j) = g(2);
    end
end

contourf(X,Y,Z,25);
% s = surf(X,Y,Z);
% s.EdgeColor = 'none';
hold on;
load('saveTable');
for i = 1:length(saveTable)
    scatter3(saveTable(i,1),...
        saveTable(i,2),...
        saveTable(i,3),...
        'MarkerEdgeColor','k','MarkerFaceColor',[1 .2 .2])
end
xlim([-extremes extremes]);
ylim([-extremes extremes]);
box on;

figure;
line = -extremes:dx:extremes;
plot(line,Z(:,floor(length(line)/2)+1));
grid minor;

for i = 1:2*extremes/dx+1
    gradNorm(i) = norm([G1(i,floor(length(line)/2)+1),...
        G2(i,floor(length(line)/2)+1)]);
end

figure;
plot(line,gradNorm,'*');
grid minor;