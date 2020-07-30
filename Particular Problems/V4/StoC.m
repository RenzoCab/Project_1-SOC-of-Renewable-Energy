function [t,v] = StoC(ulVol,olVol,delta,v0,T,V)
    t = T;
    v = (V-ulVol)*(2*delta*(1-T)+(olVol-ulVol)*T)/(olVol-ulVol)+(v0-delta)+(ulVol-(v0-delta))*T;
end