close all
clear all
clc

%% 3D:

figure;

X = [0 0 0;0 0 1;0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
v = [1;2;3;4;10;1;2;3];

[xx,yy,zz] = meshgrid(0:0.025:1);
xq = [xx(:) yy(:) zz(:)];
vq = griddatan(X,v,xq,'linear',{'QJ'});

vq = reshape(vq,size(xx));
plot3(X(:,1),X(:,2),X(:,3),'r*');
hold on;
h = slice(xx,yy,zz,vq,0:0.01:1,0:0.01:1,0:0.01:1);

set(h,'EdgeColor','none');
colorbar;
view(45,-45);
title('Interpolation (linear approximation)');

figure;

vq = griddatan(X,v,xq,'nearest',{'QJ'});

vq = reshape(vq,size(xx));
plot3(X(:,1),X(:,2),X(:,3),'r*');
hold on;
hh = slice(xx,yy,zz,vq,0:0.01:1,0:0.01:1,0:0.01:1);

set(hh,'EdgeColor','none');
colorbar;
view(45,-45);
title('Interpolation (nearest approximation)');

%% 2D:

figure;

X = [0 0; 0 1; 1 0; 1 1];
v = [1;2;3;1];

[xx,yy] = meshgrid(0:0.025:1);
xq = [xx(:) yy(:)];
vq = griddatan(X,v,xq,'linear',{'QJ'});

vq = reshape(vq,size(xx));
surf(xx,yy,vq);

set(h,'EdgeColor','none');
colorbar;
view(-65,65);
title('Interpolation (linear approximation)');

figure;

vq = griddatan(X,v,xq,'nearest',{'QJ'});

vq = reshape(vq,size(xx));
surf(xx,yy,vq);

set(h,'EdgeColor','none');
colorbar;
view(-65,65);
title('Interpolation (nearest approximation)');