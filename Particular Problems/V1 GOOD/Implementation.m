%% Constant_Demand_Compressed:

close all
clear all
clc

load('D01012017.mat');

Name_Run = 'Constant_Demand_Compressed';

MaxFlow = 4200;
MaxOil = 2000;

Vmax = 3.5005*10^9%5.0*10^9;
Vmin = 3.4995*10^9%3.0*10^9;
VminCon = 3.5*10^9%3.47*10^9;
dV = 1*10^4;
V = Vmin:dV:Vmax;
Tmax = 24;
Tmin = 1;
dt = 1;
time = Tmin:dt:Tmax;

D = 1500 + 0*500*sin(2*pi*time/24);
% D = D01012017;

u = zeros(length(V),Tmax);
u2v = zeros(length(V),Tmax);
H2s = zeros(length(V),Tmax);
TF = zeros(length(V),Tmax);
FU = zeros(length(V),Tmax);
uv = zeros(length(V),1);

eta = 8.72;
Kh = 0.0003; % Is 0.3.
Kf = 13.0; % Is 130. Here I changes 130 to 13.0 do see the effect.
H = @(V) (-7.91*10^(-7))*(V/1000000)^2 + (1.03*10^(-2))*(V/1000000) + 3.79;

A = [1 0; 0 1];
b = [MaxFlow; MaxOil];
Aeq = @(V) [eta*(H(V)-5.1)/1000 1];
beq = @(t) D(t);
lb = [0; 0];
ub = [MaxFlow; MaxOil];
f = @(i,V) [-uv(i)+Kh*eta*(H(V)-5.1) Kf];
% IT = 2.8*10^7;
IT = 2800;
op_con = zeros(length(V),2);

for t=Tmax-1:-1:1
    
    aux_Cont = op_con;
    
    for i=1:length(V)
        
        if (i~=1) && (i~=length(V)) && t < 23
            if (IT-op_con(i,1)) >= 0
                uv(i) = (u(i+1,t+1)-u(i,t+1))/(dV);
            else
                uv(i) = (u(i,t+1)-u(i-1,t+1))/(dV);
            end
        elseif i==1 && t < 23
            uv(i)=uv(i+1);
        elseif i==length(V) && t < 23
            uv(i)=uv(i-1);
        end
        
        if (V(i)>VminCon)
            b = [MaxFlow; MaxOil];
            ub = [MaxFlow; MaxOil];
        else
            b = [0; MaxOil];
            ub = [0; MaxOil];
        end
        
        op_con_out  = linprog(f(i,V(i)),A,b,Aeq(V(i)),beq(t),lb,ub);
        Ha(i) = aux_Cont(i,1)*(-uv(i)+Kh*eta*(H(V(i))-5.1)) + aux_Cont(i,2)*Kf + IT*uv(i);
        op_con(i,:) = op_con_out;
        
        TF(i,t) = op_con(i,1);
        FU(i,t) = op_con(i,2);
        
    end
    
    u2v(:,t) = uv(:);
    H2s(:,t) = Ha(:);
    u(:,t) = u(:,t+1) + Ha';
    
end

H2s(:,end) = H2s(:,end-2);
H2s(:,end-1) = H2s(:,end-2);
TF(:,end) = TF(:,end-1);
FU(:,end) = FU(:,end-1);

[v,ti] = meshgrid(V,time);
s = surf(v,ti,u');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('USD');
title('Optimal Cost Function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_1'],'epsc');

figure;
s = surf(v,ti,u2v');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('$USD/m^3$','Interpreter','latex');
title('Partial Derivative');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_2'],'epsc');

figure;
s = surf(v,ti,H2s');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('$USD/h$','Interpreter','latex');
title('Hamiltonian');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_3'],'epsc');

figure;
s = surf(v,ti,TF');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('Turbine Flow ($m^3/s$)','Interpreter','latex');
title('Turbine Flow Control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_4'],'epsc');

figure;
s = surf(v,ti,FU');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('Fuel Control ($MW/s$)','Interpreter','latex');
title('Fossil Fuel Station Control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_5'],'epsc');

figure;
p = plot(time,D);
p.LineWidth = 1;
hold on;
xlim([1 24]);
ylim([0 2000]);
grid on;
title('Constant Demand');
xlabel('Time (h)');
ylabel('Demand (MW)');
saveas(gcf,[Name_Run,'_6'],'epsc');

%% Constant_Demand:

close all
clear all
clc

load('D01012017.mat');

Name_Run = 'Constant_Demand';

MaxFlow = 4200;
MaxOil = 2000;

Vmax = 4*10^9%5.0*10^9;
Vmin = 3*10^9%3.0*10^9;
VminCon = 3.5*10^9%3.47*10^9;
dV = 1*10^7;
V = Vmin:dV:Vmax;
Tmax = 24;
Tmin = 1;
dt = 1;
time = Tmin:dt:Tmax;

D = 1500 + 0*500*sin(2*pi*time/24);
% D = D01012017;

u = zeros(length(V),Tmax);
u2v = zeros(length(V),Tmax);
H2s = zeros(length(V),Tmax);
TF = zeros(length(V),Tmax);
FU = zeros(length(V),Tmax);
uv = zeros(length(V),1);

eta = 8.72;
Kh = 0.0003; % Is 0.3.
Kf = 13.0; % Is 130. Here I changes 130 to 13.0 do see the effect.
H = @(V) (-7.91*10^(-7))*(V/1000000)^2 + (1.03*10^(-2))*(V/1000000) + 3.79;

A = [1 0; 0 1];
b = [MaxFlow; MaxOil];
Aeq = @(V) [eta*(H(V)-5.1)/1000 1];
beq = @(t) D(t);
lb = [0; 0];
ub = [MaxFlow; MaxOil];
f = @(i,V) [-uv(i)+Kh*eta*(H(V)-5.1) Kf];
% IT = 2.8*10^7;
IT = 2800;
op_con = zeros(length(V),2);

for t=Tmax-1:-1:1
    
    aux_Cont = op_con;
    
    for i=1:length(V)
        
        if (i~=1) && (i~=length(V)) && t < 23
            if (IT-op_con(i,1)) >= 0
                uv(i) = (u(i+1,t+1)-u(i,t+1))/(dV);
            else
                uv(i) = (u(i,t+1)-u(i-1,t+1))/(dV);
            end
        elseif i==1 && t < 23
            uv(i)=uv(i+1);
        elseif i==length(V) && t < 23
            uv(i)=uv(i-1);
        end
        
        if (V(i)>VminCon)
            b = [MaxFlow; MaxOil];
            ub = [MaxFlow; MaxOil];
        else
            b = [0; MaxOil];
            ub = [0; MaxOil];
        end
        
        op_con_out  = linprog(f(i,V(i)),A,b,Aeq(V(i)),beq(t),lb,ub);
        Ha(i) = aux_Cont(i,1)*(-uv(i)+Kh*eta*(H(V(i))-5.1)) + aux_Cont(i,2)*Kf + IT*uv(i);
        op_con(i,:) = op_con_out;
        
        TF(i,t) = op_con(i,1);
        FU(i,t) = op_con(i,2);
        
    end
    
    u2v(:,t) = uv(:);
    H2s(:,t) = Ha(:);
    u(:,t) = u(:,t+1) + Ha';
    
end

H2s(:,end) = H2s(:,end-2);
H2s(:,end-1) = H2s(:,end-2);
TF(:,end) = TF(:,end-1);
FU(:,end) = FU(:,end-1);

[v,ti] = meshgrid(V,time);
s = surf(v,ti,u');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('USD');
title('Optimal Cost Function');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_1'],'epsc');

figure;
s = surf(v,ti,u2v');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('$USD/m^3$','Interpreter','latex');
title('Partial Derivative');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_2'],'epsc');

figure;
s = surf(v,ti,H2s');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('$USD/h$','Interpreter','latex');
title('Hamiltonian');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_3'],'epsc');

figure;
s = surf(v,ti,TF');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('Turbine Flow ($m^3/s$)','Interpreter','latex');
title('Turbine Flow Control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_4'],'epsc');

figure;
s = surf(v,ti,FU');
xlabel('Volume ($m^3$)','Interpreter','latex');
ylabel('Time (h)');
zlabel('Fuel Control ($MW/s$)','Interpreter','latex');
title('Fossil Fuel Station Control');
xlim([Vmin Vmax])
ylim([Tmin Tmax])
box;
view(150,20);
saveas(gcf,[Name_Run,'_5'],'epsc');

figure;
p = plot(time,D);
p.LineWidth = 1;
hold on;
xlim([1 24]);
ylim([0 2000]);
grid on;
title('Constant Demand');
xlabel('Time (h)');
ylabel('Demand (MW)');
saveas(gcf,[Name_Run,'_6'],'epsc');