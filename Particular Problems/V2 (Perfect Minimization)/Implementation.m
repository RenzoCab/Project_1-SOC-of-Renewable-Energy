close all
clear all
clc

load('D01012017.mat');

MaxFlow = 4200; % In m^3/s.
MaxOil = 2000000; % In hW.

Vmax = 5.0*10^9;
Vmin = 3.0*10^9;
VminCon = 3.47*10^9;
dV = 1*10^7;
V = Vmin:dV:Vmax;
Tmax = 24*3600;
Tmin = 0;
dt = 1200;
time = Tmin:dt:Tmax;

D = 1000000*(1+sin(2*pi*time/(12*3600))); % In kW.
IT = @(t) 2800; % In m^3/s.
% D = D01012017*1000;

u = zeros(length(time),length(V));
Ph = zeros(length(time),length(V));
Pf = zeros(length(time),length(V));
PhiT = zeros(length(time),length(V));
PhiF = zeros(length(time),length(V));
uv = zeros(length(V));

eta = 8.72;
h0 = 5.1;
Khe = 0.3; % In U$S/MWh.
Kf = 130; % In U$S/MWh.
H = @(V) (-7.91*10^(-19))*(V)^2 + (1.03*10^(-8))*(V) + 3.79;
cdh = 1.37*10^(-3);
dh = @(tur) cdh*tur;
Kh = Khe*eta*(H(Vmax)-dh(MaxFlow)-5.1)/(3.6*10^6);
KH(1) = Kh;
hatKf = Kf/(3.6*10^6);
hatKh = @(uv) -uv+Kh;
K1 = @(V) eta*(H(V)-h0);
K2 = -eta*cdh;
olK2 = -K2;
Ph_star = @(V,uv) (1/(4*olK2))*(K1(V)^2-(hatKh(uv)/hatKf)^2);
olPh = @(V) eta*MaxFlow*(H(V)-dh(MaxFlow)-h0);
phit = @(V,Ph) (K1(V)-sqrt(K1(V)^2-4*olK2*Ph))/(2*olK2);
Ham = @(uv,phit,phif,t) hatKh(uv)*phit + hatKf*phif + IT(t)*uv;

uv(1)=0; % BC.
uv(length(V))=0; % BC.

for t=length(time)-1:-1:1

    for i=1:length(V)
        if (i~=1) && (i~=length(V))
            uv(i) = (u(t+1,i+1)-u(t+1,i-1))/(2*dV);
        end
        
        if (D(t)<=Ph_star(V(i),uv(i))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = D(t);
            Pf(t,i) = 0;
        elseif (D(t)<=Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i));
            if Ph_star(V(i),uv(i))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i))<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i));
            if Ph_star(V(i),uv(i))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i))>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        end
        
        PhiT(t,i) = phit(V(i),Ph(t,i));
        PhiF(t,i) = Pf(t,i);
        HAM(i) = Ham(uv(i),PhiT(t,i),PhiF(t,i),time(t));
    end
    
    u(t,:) = u(t+1,:) - HAM;

end

[v,ti] = meshgrid(V,time);

s0 = surf(v,ti,u);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('U$S');
title('Optimal cost function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(50,50);
saveas(gcf,'Result_1','epsc');

figure
grad = gradient(u);
s1 = surf(v,ti,grad);xlabel('V (m^3)');
ylabel('t (s)');
% zlabel('U$S');
title('Gradient function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(50,50);
saveas(gcf,'DU_1','epsc');

figure
s2 = surf(v,ti,PhiT);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_T');
title('Turbined flow');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,'PhiT_1','epsc');

figure
s3 = surf(v,ti,PhiF);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_F');
title('Fuel control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,'PhiF_1','epsc');

figure;
plot(time,D);
xlim([Tmin Tmax]);
grid on;
title('Demand');
xlabel('t (s)');
ylabel('D (kW)');
saveas(gcf,'D_1','epsc');

figure;
plot(time,D);
hold on;
plot(time,Ph(:,end));
plot(time,Pf(:,end));
xlim([Tmin Tmax]);
grid on;
title('Demand for V=5000 hm^3 fixed in all times');
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(size(time))*eta*4200*1000*(H(5*10^9)-h0-dh(4200))/1000,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D5000_1','epsc');

figure;
plot(time,D);
hold on;
plot(time,Ph(:,1));
plot(time,Pf(:,1));
xlim([Tmin Tmax]);
grid on;
title('Demand for V=3000 hm^3 fixed in all times');
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(size(time))*eta*4200*1000*(H(3*10^9)-h0-dh(4200))/1000,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D3000_1','epsc');

%% Otro precio de agua:

close all
clearvars -except KH
clc

load('D01012017.mat');

MaxFlow = 4200; % In m^3/s.
MaxOil = 2000000; % In hW.

Vmax = 5.0*10^9;
Vmin = 3.0*10^9;
VminCon = 3.47*10^9;
dV = 1*10^7;
V = Vmin:dV:Vmax;
Tmax = 24*3600;
Tmin = 0;
dt = 1200;
time = Tmin:dt:Tmax;

D = 1000000*(1+sin(2*pi*time/(12*3600))); % In kW. <----------XX
% D = (1+time)./(1+time)*1000000; % <----------OO 
IT = @(t) 2800; % In m^3/s.
% D = D01012017*1000;

u = zeros(length(time),length(V));
Ph = zeros(length(time),length(V));
Pf = zeros(length(time),length(V));
PhiT = zeros(length(time),length(V));
PhiF = zeros(length(time),length(V));
uv = zeros(length(V));

eta = 8.72;
h0 = 5.1;
Khe = 130; % In U$S/MWh. <----------XX
% Khe = 100; % In U$S/MWh. <----------OO
Kf = 130; % In U$S/MWh.
H = @(V) (-7.91*10^(-19))*(V)^2 + (1.03*10^(-8))*(V) + 3.79;
cdh = 1.37*10^(-3);
dh = @(tur) cdh*tur;
Kh = Khe*eta*(H(Vmax)-dh(MaxFlow)-5.1)/(3.6*10^6);
KH(2) = Kh;
hatKf = Kf/(3.6*10^6);
hatKh = @(uv) -uv+Kh;
K1 = @(V) eta*(H(V)-h0);
K2 = -eta*cdh;
olK2 = -K2;
Ph_star = @(V,uv) (1/(4*olK2))*(K1(V)^2-(hatKh(uv)/hatKf)^2);
olPh = @(V) eta*MaxFlow*(H(V)-dh(MaxFlow)-h0);
phit = @(V,Ph) (K1(V)-sqrt(K1(V)^2-4*olK2*Ph))/(2*olK2);
Ham = @(uv,phit,phif,t) hatKh(uv)*phit + hatKf*phif + IT(t)*uv;

uv(1)=0; % BC.
uv(length(V))=0; % BC.

for t=length(time)-1:-1:1

    for i=1:length(V)
        if (i~=1) && (i~=length(V))
            uv(i) = (u(t+1,i+1)-u(t+1,i-1))/(2*dV);
        end
        
        if (D(t)<=Ph_star(V(i),uv(i))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = D(t);
            Pf(t,i) = 0;
        elseif (D(t)<=Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i));
            if Ph_star(V(i),uv(i))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i))<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i));
            if Ph_star(V(i),uv(i))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i))>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        end
        
        PhiT(t,i) = phit(V(i),Ph(t,i));
        PhiF(t,i) = Pf(t,i);
        HAM(i) = Ham(uv(i),PhiT(t,i),PhiF(t,i),time(t));
    end
    
    u(t,:) = u(t+1,:) - HAM;

end

[v,ti] = meshgrid(V,time);

s0 = surf(v,ti,u);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('U$S');
title('Optimal cost function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(50,50);
saveas(gcf,'Result_2','epsc');

figure
grad = gradient(u);
s1 = surf(v,ti,grad);
xlabel('V (m^3)');
ylabel('t (s)');
% zlabel('U$S');
title('Gradient function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(-130,35);
saveas(gcf,'DU_2','epsc');

figure
s2 = surf(v,ti,PhiT,'FaceAlpha',0.5);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_T');
title('Turbined flow');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(-45,25);
saveas(gcf,'PhiT_2','epsc');

figure
s3 = surf(v,ti,PhiF);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_F');
title('Fuel control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,'PhiF_2','epsc');

figure;
plot(time,D);
xlim([Tmin Tmax]);
grid on;
title('Demand');
xlabel('t (s)');
ylabel('D (kW)');
saveas(gcf,'D_2','epsc');

figure;
plot(time,D);
hold on;
plot(time,Ph(:,end));
plot(time,Pf(:,end));
xlim([Tmin Tmax]);
grid on;
title('Demand for V=5000 hm^3 fixed in all times');
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(size(time))*eta*4200*1000*(H(5*10^9)-h0-dh(4200))/1000,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D5000_2','epsc');

figure;
plot(time,D);
hold on;
plot(time,Ph(:,1));
plot(time,Pf(:,1));
xlim([Tmin Tmax]);
grid on;
title('Demand for V=3000 hm^3 fixed in all times');
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(size(time))*eta*4200*1000*(H(3*10^9)-h0-dh(4200))/1000,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D3000_2','epsc');

%% Otro precio de agua:

close all
clearvars -except KH
clc

load('D01012017.mat');

MaxFlow = 4200; % In m^3/s.
MaxOil = 2000000; % In hW.

Vmax = 5.0*10^9;
Vmin = 3.0*10^9;
VminCon = 3.47*10^9;
dV = 1*10^7;
V = Vmin:dV:Vmax;
Tmax = 24*3600;
Tmin = 0;
dt = 1200;
time = Tmin:dt:Tmax;

% D = 1000000*(1+sin(2*pi*time/(12*3600))); % In kW. <----------XX
D = (1+time)./(1+time)*1000000; % <----------OO 
IT = @(t) 2800; % In m^3/s.
% D = D01012017*1000;

u = zeros(length(time),length(V));
Ph = zeros(length(time),length(V));
Pf = zeros(length(time),length(V));
PhiT = zeros(length(time),length(V));
PhiF = zeros(length(time),length(V));
uv = zeros(length(V));

eta = 8.72;
h0 = 5.1;
% Khe = 130; % In U$S/MWh. <----------XX
Khe = 100; % In U$S/MWh. <----------OO
Kf = 130; % In U$S/MWh.
H = @(V) (-7.91*10^(-19))*(V)^2 + (1.03*10^(-8))*(V) + 3.79;
cdh = 1.37*10^(-3);
dh = @(tur) cdh*tur;
Kh = Khe*eta*(H(Vmax)-dh(MaxFlow)-5.1)/(3.6*10^6);
KH(3) = Kh;
hatKf = Kf/(3.6*10^6);
hatKh = @(uv) -uv+Kh;
K1 = @(V) eta*(H(V)-h0);
K2 = -eta*cdh;
olK2 = -K2;
Ph_star = @(V,uv) (1/(4*olK2))*(K1(V)^2-(hatKh(uv)/hatKf)^2);
olPh = @(V) eta*MaxFlow*(H(V)-dh(MaxFlow)-h0);
phit = @(V,Ph) (K1(V)-sqrt(K1(V)^2-4*olK2*Ph))/(2*olK2);
Ham = @(uv,phit,phif,t) hatKh(uv)*phit + hatKf*phif + IT(t)*uv;

uv(1)=0; % BC.
uv(length(V))=0; % BC.

for t=length(time)-1:-1:1

    for i=1:length(V)
        if (i~=1) && (i~=length(V))
            uv(i) = (u(t+1,i+1)-u(t+1,i-1))/(2*dV);
        end
        
        if (D(t)<=Ph_star(V(i),uv(i))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = D(t);
            Pf(t,i) = 0;
        elseif (D(t)<=Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i));
            if Ph_star(V(i),uv(i))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i))<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i));
            if Ph_star(V(i),uv(i))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i))>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        end
        
        PhiT(t,i) = phit(V(i),Ph(t,i));
        PhiF(t,i) = Pf(t,i);
        HAM(i) = Ham(uv(i),PhiT(t,i),PhiF(t,i),time(t));
    end
    
    u(t,:) = u(t+1,:) - HAM;

end

[v,ti] = meshgrid(V,time);

s0 = surf(v,ti,u);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('U$S');
title('Optimal cost function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(50,50);
saveas(gcf,'Result_3','epsc');

figure
grad = gradient(u);
s1 = surf(v,ti,grad);
xlabel('V (m^3)');
ylabel('t (s)');
% zlabel('U$S');
title('Gradient function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(50,50);
saveas(gcf,'DU_3','epsc');

figure
s2 = surf(v,ti,PhiT,'FaceAlpha',0.5);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_T');
title('Turbined flow');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(-45,25);
saveas(gcf,'PhiT_3','epsc');

figure
s3 = surf(v,ti,PhiF);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_F');
title('Fuel control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,'PhiF_3','epsc');

figure;
plot(time,D);
xlim([Tmin Tmax]);
grid on;
title('Demand');
xlabel('t (s)');
ylabel('D (kW)');
saveas(gcf,'D_3','epsc');

figure;
plot(time,D);
hold on;
plot(time,Ph(:,end));
plot(time,Pf(:,end));
xlim([Tmin Tmax]);
grid on;
title('Demand for V=5000 hm^3 fixed in all times');
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(size(time))*eta*4200*1000*(H(5*10^9)-h0-dh(4200))/1000,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D5000_3','epsc');

figure;
plot(time,D);
hold on;
plot(time,Ph(:,1));
plot(time,Pf(:,1));
xlim([Tmin Tmax]);
grid on;
title('Demand for V=3000 hm^3 fixed in all times');
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(size(time))*eta*4200*1000*(H(3*10^9)-h0-dh(4200))/1000,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D3000_3','epsc');