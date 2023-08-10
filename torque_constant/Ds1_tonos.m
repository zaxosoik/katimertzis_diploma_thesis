function Ds1_tonos = Ds1_tonos(y,r,l)

% eksisosi Ds'
% Ds_tonos = ((r^2)*cos(th)*((l/r)^2-(sin(th))^2)^(3/2) + (r^2)*(sin(th))^4 - (l^2)*(sin(th))^2 + (l^2)*(cos(th))^2) / (r*((l/r)^2-(sin(th))^2)^(3/2));
% to Ds_tonos eiserxetai mono sto Tinertia

Dss_tonos = zeros(12,1);  % arxikopoiisi

% diafora fasis kilindron
n = 5;
aa = 360/n;
a = aa*pi/180;

Dss_tonos(3) = ((r^2)*cos(y(3)-pi)*((l/r)^2-(sin(y(3)))^2)^(3/2) + (r^2)*(sin(y(3)))^4 - (l^2)*(sin(y(3)))^2 + (l^2)*(cos(y(3)-pi))^2) / (r*((l/r)^2-(sin(y(3)))^2)^(3/2));
Dss_tonos(6) = ((r^2)*cos(y(6)+1*a-pi)*((l/r)^2-(sin(y(6)+1*a))^2)^(3/2) + (r^2)*(sin(y(6)+1*a))^4 - (l^2)*(sin(y(6)+1*a))^2 + (l^2)*(cos(y(6)+1*a-pi))^2) / (r*((l/r)^2-(sin(y(6)+1*a))^2)^(3/2));
Dss_tonos(5) = ((r^2)*cos(y(5)+2*a-pi)*((l/r)^2-(sin(y(5)+2*a))^2)^(3/2) + (r^2)*(sin(y(5)+2*a))^4 - (l^2)*(sin(y(5)+2*a))^2 + (l^2)*(cos(y(5)+2*a-pi))^2) / (r*((l/r)^2-(sin(y(5)+2*a))^2)^(3/2));
Dss_tonos(4) = ((r^2)*cos(y(4)+3*a-pi)*((l/r)^2-(sin(y(4)+3*a))^2)^(3/2) + (r^2)*(sin(y(4)+3*a))^4 - (l^2)*(sin(y(4)+3*a))^2 + (l^2)*(cos(y(4)+3*a-pi))^2) / (r*((l/r)^2-(sin(y(4)+3*a))^2)^(3/2));
Dss_tonos(7) = ((r^2)*cos(y(7)+4*a-pi)*((l/r)^2-(sin(y(7)+4*a))^2)^(3/2) + (r^2)*(sin(y(7)+4*a))^4 - (l^2)*(sin(y(7)+4*a))^2 + (l^2)*(cos(y(7)+4*a-pi))^2) / (r*((l/r)^2-(sin(y(7)+4*a))^2)^(3/2));

Ds1_tonos = diag(Dss_tonos);  % Ds_tonos diagonios pinakas NxN, mono ta stoixeia ton kilindron exoun ton oro Ds_tonos

% sto cos(th) thelei "- pi" gia tin sosta metafora tis gonias apo -180,180 se 0,360
