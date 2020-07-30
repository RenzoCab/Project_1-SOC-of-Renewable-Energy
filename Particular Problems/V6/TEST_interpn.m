close all;
clear all;
clc;

%% 1D:

x = [1 2 3 4];
v = [12 16 31 10];
xq = 2.5;
vq = interpn(x,v,xq,'makima');
xq2 = 1:0.1:4;
vq2 = interpn(x,v,xq2,'makima');

figure;
plot(x,v,'*',xq2,vq2,'-',xq,vq);
hold on;
plot(xq,vq,'.r','markersize',30);
legend('Samples','Cubic Interpolation','Interpolated Point');

%% 2D:

[x,y] = ndgrid((0:1));
point = [0.2,0.8];
[x1,y1] = ndgrid([0,point(1),1],[0,point(2),1]);
ff = @(x,y) (x==1) + (y==1) - 2*(x+y==2);
v = ff(x,y);
Vq = interpn(x,y,v,x1,y1,'makima');
[x2,y2] = ndgrid(0:0.01:1);
Vq2 = interpn(x,y,v,x2,y2,'makima');
figure;
s = surf(x2,y2,Vq2);
s.EdgeColor = 'none';
hold on; colorbar;
plot3(point(1),point(2),Vq(2,2),'.r','markersize',30)% ACTIVATE TO SEE THE
% DOT.

%% 3D:

[x,y,z] = ndgrid((0:1));
[x1,y1,z1] = meshgrid((0:0.005:1));
F = [0,2,3,10,5,6,0,10,3]/10;
ff = @(x,y,z) (x==0 & y==0 & z==0)*F(1) +...
    (x==0 & y==0 & z==0)*F(2) +...
    (x==1 & y==0 & z==0)*F(3) +...
    (x==0 & y==1 & z==0)*F(4) +...
    (x==1 & y==1 & z==0)*F(5) +...
    (x==0 & y==0 & z==1)*F(6) +...
    (x==1 & y==0 & z==1)*F(7) +...
    (x==0 & y==1 & z==1)*F(8) +...
    (x==1 & y==1 & z==1)*F(9);
v = ff(x,y,z);
Vq = interpn(x,y,z,v,x1,y1,z1,'linear');

slice(x1,y1,z1,Vq,0:0.5:1,0:0.5:1,0:0.5:1);
shading interp; colorbar;

%% 4D:

[w,x,y,z] = ndgrid((0:1));
p = 0.995;
point = [p,p,p,p];
[w1,x1,y1,z1] = ndgrid([0,point(1),1],[0,point(2),1],[0,point(3),1],[0,point(4),1]);
F = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
ff = @(w,x,y,z) (w==0 & x==0 & y==0 & z==0)*F(1) +...
    (w==1 & x==0 & y==0 & z==0)*F(2) +...
    (w==0 & x==1 & y==0 & z==0)*F(3) +...
    (w==1 & x==1 & y==0 & z==0)*F(4) +...
    (w==0 & x==0 & y==1 & z==0)*F(5) +...
    (w==1 & x==0 & y==1 & z==0)*F(6) +...
    (w==0 & x==1 & y==1 & z==0)*F(7) +...
    (w==1 & x==1 & y==1 & z==0)*F(8) +...
    (w==0 & x==0 & y==0 & z==1)*F(9) +...
    (w==1 & x==0 & y==0 & z==1)*F(10) +...
    (w==0 & x==1 & y==0 & z==1)*F(11) +...
    (w==1 & x==1 & y==0 & z==1)*F(12) +...
    (w==0 & x==0 & y==1 & z==1)*F(13) +...
    (w==1 & x==0 & y==1 & z==1)*F(14) +...
    (w==0 & x==1 & y==1 & z==1)*F(15) +...
    (w==1 & x==1 & y==1 & z==1)*F(16);
v = ff(w,x,y,z);
Vq = interpn(w,x,y,z,v,w1,x1,y1,z1,'makima');

Vq(2,2,2,2)
fprintf('Completed: %.2f.\n',Vq(2,2,2,2));