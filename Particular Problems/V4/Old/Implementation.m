%% The Square:

close all
clear all
clc

set(groot,'defaultFigureVisible','off') % Show plots OFF.

load('D01012017.mat');

FlowMax = 4.2*10^3; % In m^3/s.
FuelMax = 2*10^6; % In hW.
VolMax = 5*10^9; % In m^3.
VolMinCond = 3.47*10^9; % In m^3.
TMax = 24*3600; % In s.
DMax = 2*10^6; % In hW.

Vini = 0.7; % Initial condition between 0.6 and 1.
Vol = @(i) i*VolMax;

Tmax = 1;
Tmin = 0;
dt = 0.0001;
time = Tmin:dt:Tmax;

IT = @(t) 2800; % In m^3/s.
Vmax = (Vini*VolMax+IT(1)*3600*24)/VolMax;
Vmin = (Vini*VolMax+(IT(1)-FlowMax)*3600*24)/VolMax;
DV = Vmax-Vmin;
Num = 20;
dV = DV/Num;
V = Vmin:dV:Vmax;

D = DMax*(2/3+(1/3)*sin(4*pi*time)); % In kW.
% D = DMax*0.3*time./time;
IT = @(t) 2800; % In m^3/s.
% D = D01012017*1000;

u = zeros(length(time),length(V));
Ph = zeros(length(time),length(V));
Pf = zeros(length(time),length(V));
PhiT = zeros(length(time),length(V));
PhiF = zeros(length(time),length(V));
uv = zeros(1,length(V));

eta = 8.72;
h0 = 5.1;
Khe = 0.3; % In U$S/MWh.
Kf = 130; % In U$S/MWh.
H = @(V) (-19.8)*(V).^2 + (51.5)*(V) + 3.79;
cd = 5.3;
d = @(tur) cd*tur;
Kh = Khe*eta*(H(Vini)-d(1)-h0)/(3.6*10^6);
KH(1) = Kh;
hatKf = Kf*FuelMax/(3.6*10^6);
hatKh = @(uv) FlowMax*(Kh-uv/VolMax);
K1 = @(V) eta*FlowMax*(H(V)-h0);
K2 = -eta*cd*FlowMax;
olK2 = -K2;
Ph_star = @(V,uv) (1/(4*olK2))*(K1(V).^2-(hatKh(uv)*FuelMax/hatKf).^2);
olPh = @(V) eta*FlowMax*(H(V)-d(1)-h0);
phit = @(V,Ph) (K1(V)-sqrt(K1(V)^2-4*olK2*Ph))/(2*olK2);
Ham = @(uv,phit,phif,t) TMax*(hatKh(uv)*phit + hatKf*phif + IT(t)*uv/VolMax);
% D = ones(1,length(time))*olPh(Vmax);

uv(1)=0; % BC.
uv(length(V))=0; % BC.

for t=length(time):-1:1

    for i=1:length(V)
        if (i~=1) && (i~=length(V)) && (t~=length(time))
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
        PhiF(t,i) = Pf(t,i)/FuelMax;
        HAM(i) = Ham(uv(i),PhiT(t,i),PhiF(t,i),time(t));
    end
    
    if t~=length(time)
        u(t,:) = u(t+1,:) + dt*HAM;
    end

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
s0.EdgeColor = 'none';
view(-225,50);
saveas(gcf,'Result_1','epsc');

figure
grad = gradient(u);
s1 = surf(v,ti,grad);
xlabel('V (m^3)');
ylabel('t (s)');
title('Gradient function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
s1.EdgeColor = 'none';
view(-225,50);
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
s2.EdgeColor = 'none';
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
s3.EdgeColor = 'none';
saveas(gcf,'PhiF_1','epsc');

figure;
plot(time,D);
xlim([Tmin Tmax]);
grid on;
title('Demand in normalized time');
xlabel('t (normalized)');
ylabel('D (kW)');
saveas(gcf,'D_1','epsc');

figure;
hold on;
grid on;
title(['Demand for V=',num2str(Vmax*VolMax/1000000),' hm^3 fixed in all times']);
xlabel('t (s)');
ylabel('D (kW)');
plot(time,D/DMax);
plot(time,Ph(:,end)/DMax);
plot(time,Pf(:,end)/DMax);
plot(time,ones(1,length(time))*olPh(Vmax)/DMax,'o');
xlim([Tmin Tmax]);
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D5000_1','epsc');

% figure;
% hold on;
% grid on;
% title('Demand for V=4000 hm^3 fixed in all times');
% xlabel('t (s)');
% ylabel('D (kW)');
% plot(time,D/DMax);
% plot(time,Ph(:,floor(end/2))/DMax);
% plot(time,Pf(:,floor(end/2))/DMax);
% plot(time,ones(1,length(time))*olPh((Vmax-Vmin)/2+Vmin)/DMax,'o');
% xlim([Tmin Tmax]);
% legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
% saveas(gcf,'D4000_1','epsc');

figure;
plot(time,D/DMax);
hold on;
plot(time,Ph(:,1)/DMax);
plot(time,Pf(:,1)/DMax);
xlim([Tmin Tmax]);
grid on;
title(['Demand for V=',num2str(Vmin*VolMax/1000000),' hm^3 fixed in all times']);
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(1,length(time))*olPh(Vmin)/DMax,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
saveas(gcf,'D3000_1','epsc');

Data1 = {V,time,u};

%% The Cone:

close all
clearvars -except Data1 Vini Vmax Vmin IT
clc

FlowMax = 4.2*10^3; % In m^3/s.
FuelMax = 2*10^6; % In hW.
VolMax = 5*10^9; % In m^3.
TMax = 24*3600; % In s.
DMax = 2*10^6; % In hW.

de = 0.01; % Delta.

Vol = @(i) i*VolMax;

Tmax = 1;
Tmin = 0;
dt = 0.01;
time = Tmin:dt:Tmax;

D = DMax*(2/3+(1/3)*sin(4*pi*time)); % In kW.
DV = Vmax-Vmin;
Num = 20;
dV = DV/Num;
V = Vmin:dV:Vmax;

u = zeros(length(time),length(V));
Ph = zeros(length(time),length(V));
Pf = zeros(length(time),length(V));
PhiT = zeros(length(time),length(V));
PhiF = zeros(length(time),length(V));
uv = zeros(1,length(V));

eta = 8.72;
h0 = 5.1;
Khe = 0.3; % In U$S/MWh.
Kf = 130; % In U$S/MWh.
H = @(V) (-19.8)*(V)^2 + (51.5)*(V) + 3.79;
cd = 5.3;
d = @(tur) cd*tur;
Kh = Khe*eta*(H(Vini)-d(1)-h0)/(3.6*10^6);
KH(2) = Kh;
hatKf = Kf*FuelMax/(3.6*10^6);
K1 = @(V) eta*FlowMax*(H(V)-h0);
K2 = -eta*cd*FlowMax;
olK2 = -K2;
olPh = @(V) eta*FlowMax*(H(V)-d(1)-h0);
phit = @(V,Ph) (K1(V)-sqrt(K1(V)^2-4*olK2*Ph))/(2*olK2);
Fv = @(t) DV/(2*de*(1-t)+DV*t);
hatKh = @(uv,t) FlowMax*(Kh-uv*Fv(t)/VolMax);
Ham = @(uv,phit,phif,t) TMax*(hatKh(uv,t)*phit + hatKf*phif + IT(t)*uv*Fv(t)/VolMax);
Ph_star = @(V,uv,t) (1/(4*olK2))*(K1(V)^2-(hatKh(uv,t)*FuelMax/hatKf)^2);
Ft = @(t,v) DV*(-Vmin-de+Vini)/(2*de*(1-t)+DV*t)-(DV*(DV-2*de)*(-t*(Vmin+de-Vini)+de+v-Vini))/(2*de*(1-t)+DV*t)^2;
% D = ones(1,length(time))*olPh(Vmax);

uv(1)=0; % BC.
uv(length(V))=0; % BC.

for t=length(time):-1:1

    for i=1:length(V)
        if (i~=1) && (i~=length(V)) && (t~=length(time))
            uv(i) = (u(t+1,i+1)-u(t+1,i-1))/(2*dV);
        end
        
        if (D(t)<=Ph_star(V(i),uv(i),time(t))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = D(t);
            Pf(t,i) = 0;
        elseif (D(t)<=Ph_star(V(i),uv(i),time(t))) && (D(t)>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        elseif (D(t)>Ph_star(V(i),uv(i),time(t))) && (D(t)<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i),time(t));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i),time(t));
            if Ph_star(V(i),uv(i),time(t))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i),time(t))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i),time(t))<=olPh(V(i)))
            Ph(t,i) = Ph_star(V(i),uv(i),time(t));
            Pf(t,i) = D(t)-Ph_star(V(i),uv(i),time(t));
            if Ph_star(V(i),uv(i),time(t))<0
                Ph(t,i) = 0;
                Pf(t,i) = D(t);
            end
        elseif (D(t)>Ph_star(V(i),uv(i),time(t))) && (D(t)>olPh(V(i))) && (Ph_star(V(i),uv(i),time(t))>olPh(V(i)))
            Ph(t,i) = olPh(V(i));
            Pf(t,i) = D(t)-olPh(V(i));
        end
        PhiT(t,i) = phit(V(i),Ph(t,i));
        PhiF(t,i) = Pf(t,i)/FuelMax;
        HAMtot(i) = Ham(uv(i),PhiT(t,i),PhiF(t,i),time(t))+uv(i)*Ft(time(t),V(i));
    end
    
    if t~=length(time)
        u(t,:) = u(t+1,:) + dt*HAMtot;
    end

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
s0.EdgeColor = 'none';
% saveas(gcf,'Result_1','epsc');

figure
grad = gradient(u);
s1 = surf(v,ti,grad);
xlabel('V (m^3)');
ylabel('t (s)');
title('Gradient function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
s1.EdgeColor = 'none';
% saveas(gcf,'DU_1','epsc');

figure
s2 = surf(v,ti,PhiT);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_T');
title('Turbined flow');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
s2.EdgeColor = 'none';
% saveas(gcf,'PhiT_1','epsc');

figure
s3 = surf(v,ti,PhiF);
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('\phi_F');
title('Fuel control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
s3.EdgeColor = 'none';
% saveas(gcf,'PhiF_1','epsc');

figure;
plot(time,D);
xlim([Tmin Tmax]);
grid on;
title('Demand in normalized time');
xlabel('t (normalized)');
ylabel('D (kW)');
% saveas(gcf,'D_1','epsc');

figure;
hold on;
grid on;
title(['Demand for V=',num2str(Vmax*VolMax/1000000),' hm^3 fixed in all times']);
xlabel('t (s)');
ylabel('D (kW)');
plot(time,D/DMax);
plot(time,Ph(:,end)/DMax);
plot(time,Pf(:,end)/DMax);
plot(time,ones(1,length(time))*olPh(Vmax)/DMax,'o');
xlim([Tmin Tmax]);
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
% saveas(gcf,'D5000_1','epsc');

% figure;
% hold on;
% grid on;
% title('Demand for V=4000 hm^3 fixed in all times');
% xlabel('t (s)');
% ylabel('D (kW)');
% plot(time,D/DMax);
% plot(time,Ph(:,floor(end/2))/DMax);
% plot(time,Pf(:,floor(end/2))/DMax);
% plot(time,ones(1,length(time))*olPh((Vmax-Vmin)/2+Vmin)/DMax,'o');
% xlim([Tmin Tmax]);
% legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
% % saveas(gcf,'D4000_1','epsc');

figure;
plot(time,D/DMax);
hold on;
plot(time,Ph(:,1)/DMax);
plot(time,Pf(:,1)/DMax);
xlim([Tmin Tmax]);
grid on;
title(['Demand for V=',num2str(Vmin*VolMax/1000000),' hm^3 fixed in all times']);
xlabel('t (s)');
ylabel('D (kW)');
plot(time,ones(1,length(time))*olPh(Vmin)/DMax,'o');
legend('D(t)','P_H(t)','P_F(t)','P_{H_{MAX}}');
% saveas(gcf,'D3000_1','epsc');

Data2 = {V,time,u};

%% Comparison:

IWant = 0;
if IWant == 1

close all
clearvars -except Data1 Data2 Vini de Vmax Vmin
clc

all = 1; % If =1 we see all, if =2 we see half.

figure;
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('U$S');
title('Optimal cost function comparison');
xlim([Vmin Vmax])
ylim([0 1])
hold on;
box;
grid on;
for i=1:all:length(Data1{1})
    for j=1:all:length(Data1{2})
        
        scatter3(Data1{1}(i),Data1{2}(j),Data1{3}(j,i),'MarkerFaceColor',[.75 0 .75]);
        
    end
end

for i=1:length(Data2{2}) % Time.
    for j=1:length(Data2{1}) % Volume.
        
        [~,S2C(i,j)] = StoC(Vmin,Vmax,de,Vini,Data2{2}(i),Data2{1}(j));
        
    end
end

for i=1:all:length(Data2{1}) % Volume.
    for j=1:all:length(Data2{2}) % Time.
        
        scatter3(S2C(j,i),Data2{2}(j),Data2{3}(j,i),'MarkerFaceColor',[0 .75 .75]);
        
    end
end

% saveas(gcf,'Comparasion1','epsc');
end

%% Comparison (another plot):

set(groot,'defaultFigureVisible','on') % Show plots ON.
close all
clearvars -except Data1 Data2 Vini de Vmax Vmin
clc

all = 1; % If =1 we see all, if =2 we see half.

figure;
[v,ti] = meshgrid(Data1{1},Data1{2});
s0 = surf(v,ti,Data1{3});
xlabel('V (m^3)');
ylabel('t (s)');
zlabel('U$S');
title('Optimal cost function comparison');
xlim([Vmin Vmax])
ylim([0 1])
s0.EdgeColor = 'none';
hold on;
box;
grid on;

for i=1:length(Data2{2}) % Time.
    for j=1:length(Data2{1}) % Volume.
        
        [~,S2C(i,j)] = StoC(Vmin,Vmax,de,Vini,Data2{2}(i),Data2{1}(j));
        
    end
end

for i=1:all:length(Data2{1}) % Volume.
    for j=1:all:length(Data2{2}) % Time.
        
        scatter3(S2C(j,i),Data2{2}(j),Data2{3}(j,i),'MarkerFaceColor',[0 .75 .75]);
        
    end
end

view(145,30);
saveas(gcf,'Comparison2','epsc');
saveas(gcf,'Comparison2.png');