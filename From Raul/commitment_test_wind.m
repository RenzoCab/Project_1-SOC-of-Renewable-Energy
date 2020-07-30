% test for the single commitment of a thermal unit.
% small test adding uncertain wind.
% for simplicity, the maximum value of the wind power is set to the minimum
% demand value.

close all, clear all
% Beware, the initial time is always zero, the final time is T
hFD = 1e-3; %used for numerical (Finite Difference) differentiation.
ihFD = 1/hFD;

CASED = 2;
switch CASED
    case 1
        T = 1; % final time
        N = 100; % Time steps
        Nw = 10; % Intervals in the wind direction
        
        % Thermal Unit definition
        Cost_Start = 100; %400; % Only for the first thermal unit
        Var_Cost = [5;15]; % first is for the commitment unit, second is for no commitment unit
        Max_Pow  = 200; % only bound for the commitment unit;
        Min_Pow = 90; %only bound for the commitment unit
        DvarPow = Max_Pow-Min_Pow;
        
        % Instantaneous Power Demand definition
        tt = linspace(0,T,N+1);
        Dmax = 400; % Max demand
        Dmin = Dmax/10;
        Demand = Dmin+(Dmax-Dmin)*(tt/T).^4;
        
        % Wind power grid definition
        
        Pwind_max = Dmin;
        Pwindg = linspace(0,1,Nw+1)'; % work in the interval [0,1], then rescale
        
        % Wind power model: assumes constant forecast
        theta = 4*T;
        alpha = 0.5;
        
        P_forecastf = inline('0.1+0.8*sin(pi*t/T).^4','t','T');

        %P_forecastf = inline('0.5*ones(size(t))','t','T');
        % time independent diffusion values, the advection is defined later
        % calling the forecast
        diffu_vals = 0.5*(alpha*theta.*Pwindg.*(1-Pwindg)); %b^2/2, it is time independent so it is computed here
        
        
    case 2
        T = 1; % final time
        N = 100; % Time steps
        Nw = 10; % Intervals in the wind direction
        
        % Thermal Unit definition
        Cost_Start = 100; %400; % Only for the first thermal unit
        Var_Cost = [5;10]; % first is for the commitment unit, second is for no commitment unit
        Max_Pow  = 200; % only bound for the commitment unit;
        Min_Pow = 100; %only bound for the commitment unit
        DvarPow = Max_Pow-Min_Pow;
        
        % Instantaneous Power Demand definition
        tt = linspace(0,T,N+1);
        Dmax = 400; % Max demand
        Dmin = Dmax/7;
        Demand = Dmin+(Dmax-Dmin)*sin(2*pi*tt/T).^4;
        
        % Wind power grid definition
        
        Pwind_max = Dmin;
        Pwindg = linspace(0,1,Nw+1)'; % work in the interval [0,1], then rescale
        
        % Wind power model: assumes constant forecast
        theta = 12*T;
        alpha = 0.4;
        P_forecastf = inline('0.2+0.6*sin(2*pi*t/T).^4','t','T');
        diffu_vals = 0.5*(alpha*theta.*Pwindg.*(1-Pwindg)); %b^2/2
        
        
    case 3
        
        T = 1; % final time
        N = 100; % Time steps
        Nw = 5; % Intervals in the wind direction
        
        % Thermal Unit definition
        Cost_Start = 100; %400; % Only for the first thermal unit
        Var_Cost = [5;10]; % first is for the commitment unit, second is for no commitment unit
        Max_Pow  = 200; % only bound for the commitment unit;
        Min_Pow = 90; %only bound for the commitment unit
        DvarPow = Max_Pow-Min_Pow;
        
        % Instantaneous Power Demand definition
        tt = linspace(0,T,N+1);
        Dmax = 400; % Max demand
        Dmin = Dmax/10;
        Demand = Dmin+(Dmax-Dmin)*(tt/T-1).^4;
        
        
        % Wind power grid definition
        
        Pwind_max = Dmin;
        Pwindg = linspace(0,1,Nw+1)'; % work in the interval [0,1], then rescale
        % Wind power model: assumes constant forecast
        theta = 12*T;
        alpha = 1;
        P_forecastf = inline('0.1+0.8*sin(pi*t/T).^4','t','T');
        diffu_vals = 0.5*(alpha*theta.*Pwindg.*(1-Pwindg)); %b^2/2
        
        
        %%%%%%% END OF CASES DEFINITION %%%%%%%
    otherwise
        error('Input case not defined yet, please verify')
        
        
        
end %switch

% function X = d_P_forecast_dtf(t,T)
% % forecasted percentage of wind power
% X = 2*ihFD*(P_forecastf(t+hFD,T)-P_forecastf(t-hFD,T));
% end %d_P_forecast_dtf


dt = T/N; % time step
U = zeros(N+1,Nw+1,2); % cost to go function, (:,1)=OFF state; (:,2) = ON state

DPow = 1/Nw; % grid size in wind direction

%d_P_forecast_vals_dt = diff(P_forecast_vals)./diff(tt); % numerical approximation of time derivative
%d_P_forecast_vals_dt = [d_P_forecast_vals_dt d_P_forecast_vals_dt(end)]; % trivial extrapolation to match sizes


%%%%%%%%%%%%%%%%%%%%%%%%%

% diffusion operator*b^2/2, centered differences. No need for b.c., the
% diffusion is zero at the boundaries.
e = diffu_vals/(DPow^2);
Diffu_op = spdiags(-2*e, 0, Nw+1, Nw+1);
Diffu_op = spdiags([0;e(1:Nw)],1,Diffu_op);
Diffu_op = spdiags(e(2:Nw+1),-1,Diffu_op);



%keyboard
for mm=N:-1:1 % Dynamic programming Backward iteration, split method with the wind power generator
    
    % advection part
t = tt(mm);
drift_vals = 0.5*ihFD*(P_forecastf(t+hFD,T)-P_forecastf(t-hFD,T))-theta * (Pwindg-P_forecastf(t,T)); % mean reverting drift, this is the latest model that Waleed is using.

Advec = Advec_op(drift_vals);
%
LKB2 = -Advec + Diffu_op;

LKBaux = speye(Nw+1,Nw+1)-LKB2*dt; % for later time stepping
    
    % Optimal controls
    Demand_reduced = max(Demand(mm)-Pwind_max*Pwindg',0); % Demand - Wind Power
    phis = max(min(Demand_reduced-Min_Pow,DvarPow),0); % tentative variable power of committed unit
    P2 = max(Demand_reduced-phis-Min_Pow,0); % rest of the power fulfilled by the simpler unit if commited unit is ON
    % split method, first take a minimization step, paramettric in the wind
    % step 1: backward iteration for the control
    
    Transition_Cost_to_OFF = Var_Cost(2)*Demand_reduced*dt;
    Transition_Cost_ON_to_ON = (Var_Cost(1)*(Min_Pow+phis)+Var_Cost(2)*P2)*dt;
    Transition_Cost_OFF_to_ON = Cost_Start + Transition_Cost_ON_to_ON;
    
    Future_Cost_to_OFF = Transition_Cost_to_OFF+U(mm+1,:,1);
    
    U(mm,:,1) = min([Future_Cost_to_OFF;Transition_Cost_OFF_to_ON+U(mm+1,:,2)]);
    pest = Transition_Cost_ON_to_ON +U(mm+1,:,2);
    U(mm,:,2) = min([Future_Cost_to_OFF;pest]);
    % step 2: backward iteration in the wind direction
    U(mm,:,1) = (LKBaux\U(mm,:,1)')';
    U(mm,:,2) = (LKBaux\U(mm,:,2)')';
end

figure, plot(tt(1:end-1),diff(U(:,:,1)),'k'), hold on,  plot(tt(1:end-1),diff(U(:,:,2)),'r'), xlabel('time'), ylabel('time forward difference of cost to go function'), grid,
figure, plot(tt,U(:,:,1),'k'), hold on,  plot(tt,U(:,:,2),'r'), xlabel('time'), ylabel('cost to go function'), legend('OFF state','ON state'), grid,
figure, plot(tt,U(:,:,1)-U(:,:,2),'k'),  xlabel('time'), ylabel('OFF-ON cost to go function difference'), hold on, plot(tt,Cost_Start*ones(1,N+1),'r'),  grid,

% Forward iteration to define optimal path (primal feasible)
state = zeros(N+1,1); % state trajectory
state(1) = 1; % Initial condition: 1=OFF; 2 =ON;
Running_Cost = zeros(N+1,1); % Accumulated cost until this time
Committed_Power = zeros(N,1); % Power produced by commitment unit


% Generate a forward path realization of the wind power.
XWind = ones(N+1,1)*P_forecastf(0,T); % deterministic trajectory to test first
DeltaW = randn(N,1)*sqrt(dt);
Demand_reduced = zeros(N+1,1);
Power_Excess = zeros(N,1);
Xzerodrift = ones(N,1);

for mm=1:N
    
    
    t = tt(mm);
    drift_vals = 0.5*ihFD*(P_forecastf(t+hFD,T)-P_forecastf(t-hFD,T))-theta * (Pwindg-P_forecastf(t,T)); % mean reverting drift, this is the latest model that Waleed is using.

    Advec = Advec_op(drift_vals);
    %if phys_time>0.4, disp('check'), keyboard, end
    % Optimal controls
    Demand_reduced(mm) = max(Demand(mm)-XWind(mm)*Pwind_max,0); % Demand - Wind Power
    phis = max(min(Demand_reduced(mm)-Min_Pow,DvarPow),0); % tentative variable power of committed unit
    P2 = max(Demand_reduced(mm)-phis-Min_Pow,0); % rest of the power fulfilled by the simpler unit if commited unit is ON
    Power_Excess(mm) = max(XWind(mm)*Pwind_max+Min_Pow*(state(mm)-1)-Demand(mm),0);
    
    Transition_Cost_to_OFF = Var_Cost(2)*Demand_reduced(mm)*dt;
    Transition_Cost_ON_to_ON = (Var_Cost(1)*(Min_Pow+phis)+Var_Cost(2)*P2)*dt;
    Transition_Cost_OFF_to_ON = Cost_Start+Transition_Cost_ON_to_ON;
    
    %    interpolate Delta_U(mm+1,wind,1:2) based on the generator of the wind SDE
    %    compute mathcalL(U) for both ON and OFF options, (partial_t U+ 1/2b^2Dxx U + a(x) Dx U) * dt
    Delta_U_1 = U(mm+1,:,1)-U(mm,:,1); % partial_t U * dt contribution
    Delta_U_1 = Delta_U_1 + ((-Advec+Diffu_op)*U(mm,:,1)')'*dt; % advection-diffusion component
    Delta_U_1_current = interp1(Pwindg,Delta_U_1,XWind(mm),'spline');
    
    Delta_U_2 = U(mm+1,:,2)-U(mm,:,2);
    Delta_U_2 = Delta_U_2 + ((-Advec+Diffu_op)*U(mm,:,2)')'*dt; % advection-diffusion component
    Delta_U_2_current = interp1(Pwindg,Delta_U_2,XWind(mm),'spline');
    
    if state(mm)==1  % we are currently OFF
        %
        %here is the mistake, it is in the evaluation of the alternatives
        Cost_alternative = [Transition_Cost_to_OFF+Delta_U_1_current;
            Transition_Cost_OFF_to_ON+interp1(Pwindg,U(mm,:,2)-U(mm,:,1),XWind(mm),'spline')+Delta_U_2_current];
        [aux,state(mm+1)] = min(Cost_alternative);
        if state(mm+1)==1 % we continue OFF
            Delta_Cost = Transition_Cost_to_OFF;
            Committed_Power(mm) =0;
        else % we switch to ON
            Committed_Power(mm) = Min_Pow+phis;
            Delta_Cost = Transition_Cost_OFF_to_ON;
        end
    else % state(mm)=2, % we are currently ON
        Cost_alternative = [Transition_Cost_to_OFF+interp1(Pwindg,U(mm,:,1)-U(mm,:,2),XWind(mm),'spline')+Delta_U_1_current;
            Transition_Cost_ON_to_ON+Delta_U_2_current];
        [aux,state(mm+1)] = min(Cost_alternative);
        if state(mm+1)==1 % we switch to OFF
            Committed_Power(mm) =0;
            Delta_Cost = Transition_Cost_to_OFF;
        else % we continue ON
            Committed_Power(mm) = Min_Pow+phis;
            Delta_Cost = Transition_Cost_ON_to_ON;
        end
    end
    Running_Cost(mm+1) = Running_Cost(mm) + Delta_Cost;
    
    % time stepping to simulate XWind
    current_drift = ... %interp1(Pwindg,drift_vals,XWind(mm)); % -theta * (Pwindg-P_forecast)
    0.5*ihFD*(P_forecastf(t+hFD,T)-P_forecastf(t-hFD,T))-theta * (XWind(mm)-P_forecastf(t,T));
    current_diffu = interp1(Pwindg,sqrt(2*diffu_vals),XWind(mm)); % beware of redefinition sqrt(2*(1/2*b^2))
%     if or((XWind(mm)==0),(XWind(mm)==1)) 
%         disp 'check',
%         t,
%         0.5*ihFD*(P_forecastf(t+hFD,T)-P_forecastf(t-hFD,T)),
%         -theta*(XWind(mm)-P_forecastf(t,T)),
%         figure,
%         plot(tt,P_forecastf(tt,T))
%         hold on, plot(tt,XWind) 
%         keyboard, 
%     end
    XWind(mm+1) = max(min(XWind(mm)+dt*current_drift+ current_diffu*DeltaW(mm),1),0); % projected Euler Maruyama
    Xzerodrift(mm)= 0.5*ihFD*(P_forecastf(t+hFD,T)-P_forecastf(t-hFD,T))/theta+ P_forecastf(t+hFD,T);
end
figure; plot(tt,state); xlabel 'time', ylabel 'Commitment state'
figure; plot(tt,Running_Cost); xlabel 'time', ylabel('Running Cost'),grid
figure; plot(tt(1:N),XWind(1:N)); xlabel 'time',  hold on, plot(tt,P_forecastf(tt,T));plot(tt(1:N),Xzerodrift); legend('Wind Power','Forecast','Xzero'), grid

figure; plot(tt,Demand), hold on;grid
plot(tt(1:end-1),Committed_Power,'r'),
plot(tt(1:end-1),Committed_Power+XWind(1:N)*Pwind_max,'b'), 
xlabel 'time', ylabel 'Power balance',
legend('Demand','Committed Power','Committed Power+Wind Power')
figure; plot(tt(1:N),Power_Excess),grid
xlabel 'time', ylabel 'Power Excess',

function Advec = Advec_op(drift_vals)

Nw = length(drift_vals)-1;
DPow = 1/Nw;

aux = drift_vals/DPow;
upwind = 1; %upwinding scheme for the drift term
if upwind
    LKB1 = speye(Nw+1);
    LKB1 = spdiags(-ones(Nw+1,1),1,LKB1);
    LKB1(Nw+1,Nw+1) = -1;
    LKB1(Nw+1,Nw) = 1;
    
    LKB2 = -speye(Nw+1);
    LKB2 = spdiags(ones(Nw+1,1),-1,LKB2);
    LKB2(1,2) = -1;
    LKB2(1,1) = 1;
    
    iip = find(drift_vals>0);liip = length(iip);
    iin = find(drift_vals<=0);liin = length(iin);
    
    Advec = sparse(Nw+1,Nw+1);
    Advec(iip,:) = spdiags(aux(iip),0,liip,liip)* LKB1(iip,:);
    Advec(iin,:) = spdiags(aux(iin),0,liin,liin)* LKB2(iin,:);
    
else
    % test with centered differences, periodic boundary conditions, all
    % imaginary eigenvalues
    LKB = sparse(Nw+1,Nw+1);
    LKB = spdiags(ones(Nw+1,1),1,LKB);
    LKB = spdiags(-ones(Nw+1,1),-1,LKB);
    LKB(1,Nw+1) = -1;
    LKB(Nw+1,1) = 1;
    %%%%%%
    % one sided difference, periodic boundary conditions (all eigenvalues with
    % positive real part)
    LKB = sparse(Nw+1,Nw+1);
    LKB = spdiags(ones(Nw+1,1),0,LKB);
    LKB = spdiags(-ones(Nw+1,1),-1,LKB);
    LKB(1,Nw+1) = -1;
    %%%%%%
    
    % LKB = sparse(Nw+1,Nw+1);
    % LKB = spdiags([0;aux(1:Nw)],1,LKB);
    % LKB = spdiags(-aux(2:Nw+1),-1,LKB);
    % LKB(1,1) = -aux(1);
    % LKB(Nw+1,Nw+1) = -aux(Nw+1);
end
%full(LKB)
%pause

end %function