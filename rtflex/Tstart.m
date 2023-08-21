function Tstart = Tstart(y,Ds)  % ropi logo kaysis

% ypologismos eksisosis piesis
Pstart = zeros(12,1);  % arxikopoiisi

%if y(5) < pi+37*2*pi/360
if y(5) < 115*2*pi/360 + 37*2*pi/360

    Pstart(5) = 30*10^5;   % Pascal (30 bar), statheri piesi
else
    
    Pstart(5) = 0;
end


B = 0.5;     % bore (m)
l = 2.05;   % piston rod
r = 2.05/2;  % crankshaft, stroke/2

Tstart = - Ds*(pi/4)*(B^2)*Pstart;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
end

% sxolia
% theoreitai sinexis piesi 30 bar kai askeitai gia maximum 115 moires h
% otan anoigei h exhaust valve (opoio einai mikrotero)
% opening exhaust valve 130 moires
% rpm pou pianei sto telos tou starting air system : 28.7 rpm