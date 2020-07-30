function [t,v] = CtoS(ulVol,olVol,delta,v0,T,V)
    t = T;
    v = ulVol+(olVol-ulVol)/(2*delta*(1-T)+(olVol-ulVol)*T)*(V-(v0-delta)-(ulVol-(v0-delta))*T);
end