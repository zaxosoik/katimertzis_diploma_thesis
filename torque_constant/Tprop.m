function Tprop = Tprop(wmega)  % ropi fortiou elikas

PMCR = 6350*10^3;          % Watt
NMCR = 99;                 % rpm
WMCR = 2*pi*NMCR/60;       % rad/s
Tp = PMCR*wmega^2/WMCR^3;  % eksisosi ropis nomou elikas
%Tp = 0.85*PMCR*wmega^2/(0.97*WMCR)^3;  % eksisosi ropis nomou elikas, clean hull

Tprop = [0,0,0,0,0,0,0,0,0,0,0,Tp]';   % mono sto teleutaio DoF askeitai h ropi tis elikas

end