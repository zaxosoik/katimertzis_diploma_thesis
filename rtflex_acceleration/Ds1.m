function Ds1 = Ds1(y,r,l)

% Ds = r*sin(th)*(1+cos(th)/sqrt((l/r)^2-(sin(th))^2));  eksisosi Ds

Dss = zeros(12,1);   % arxikopoiisi

% diafora fasis kilindron
n = 5;
aa = 360/n;
a = aa*pi/180;

Dss(3) = r*sin(y(3)).*(1+cos(y(3)-pi)/sqrt((l/r)^2-(sin(y(3))).^2));
Dss(6) = r*sin(y(6)+1*a).*(1+cos(y(6)+1*a-pi)/sqrt((l/r)^2-(sin(y(6)+1*a)).^2));
Dss(5) = r*sin(y(5)+2*a).*(1+cos(y(5)+2*a-pi)/sqrt((l/r)^2-(sin(y(5)+2*a)).^2));
Dss(4) = r*sin(y(4)+3*a).*(1+cos(y(4)+3*a-pi)/sqrt((l/r)^2-(sin(y(4)+3*a)).^2));
Dss(7) = r*sin(y(7)+4*a).*(1+cos(y(7)+4*a-pi)/sqrt((l/r)^2-(sin(y(7)+4*a)).^2));

Ds1 = diag(Dss);    % Ds diagonios pinakas NxN, mono ta stoixeia ton kilindron exoun ton oro Ds


% sto cos(th) thelei "- pi" gia tin sosta metafora tis gonias apo -180,180 se 0,360
