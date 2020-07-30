function [c,ceq] = nonlcon(x)

    global Demand VG1 VG4
    eta = 8.72;
    h04 = 5.1;
    h01 = 54;
    H4 = @(V) (-19.8)*(V).^2 + (51.5)*(V) + 3.79;
    H1 = @(V) (-3.74)*(V).^2 + (16.7)*(V) + 67.7;
    cd4 = 5.3;
    d4 = @(tur) cd4*tur;
    cd1 = 1.0;
    d1 = @(tur) cd1*tur;
    FlowMax4 = 4.2*10^3; % In m^3/s.
    FlowMax1 = 640; % In m^3/s.
    FuelMax = 2*10^6; % In hW.
    
    PowH1 = @(phihatH1) eta*FlowMax1*phihatH1*(H1(VG1)-d1(phihatH1)-h01);
    PowH4 = @(phihatH4) eta*FlowMax4*phihatH4*(H4(VG4)-d4(phihatH4)-h04);
    PowF = @(phihatF) FuelMax*phihatF;
    
    c = PowH1(x(1))+PowH4(x(2))+PowF(x(3))-Demand;
    ceq = PowH1(x(1))+PowH4(x(2))+PowF(x(3))-Demand;

end