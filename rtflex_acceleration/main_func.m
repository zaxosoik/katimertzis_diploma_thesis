function main_func = main_func(time,y)

% N plithos DoF, N=12

% idioi pinakes inertia, stiffness, damping kai sto main (dedomena gia 5RTflex50D)
% SOS : pinakas y  :  sthles 1:N gonies kathe DoF, sthles N+1:2*N goniakes taxitites kathe DoF
% y(1:N) gonies, y(N+1:2*N) goniakes taxitites
% oi seires tou pinaka y deixnoun tin metavoli sto xrono ton gonion kai ton goniakon taxititon


% inertia, pinakas NxN diagonios
JJ = [8730 578 4920 4920 4920 4920 4920 735 8136 82.095 194.14 21497];
J = diag(JJ);

% stiffness, pinakas NxN tridiagonios
KK1 = [18*10^6 1.508*10^9 1.344*10^9 1.344*10^9 1.344*10^9 1.344*10^9 0.892*10^9 0.509*10^9 0.02*10^9 62.058*10^9 19.377*10^9];
K2 = -diag(KK1(1:11),-1);
K3 = -diag(KK1(1:11),1);
K4 = K2 + K3;
K1 = zeros(length(KK1)+1,1);
K1(1) = KK1(1);
K1(end)= KK1(end);

for i = 2:length(KK1)
    K1(i) = abs(K2(i,i-1)+K3(i,i+1));
end
K1 = diag(K1);
    
K = K1 + K4;

% damping, pinakas NxN tridiagonios
CC1 = [200000 0 0 0 0 0 0 0 0 0 0];
C2 = -diag(CC1(1:11),-1);
C3 = -diag(CC1(1:11),1);
C4 = C2 + C3;
C1 = zeros(length(CC1)+1,1);
C1(1) = CC1(1);
C1(end) = CC1(end);

for i = 2:length(CC1)
    C1(i) = abs(C2(i,i-1)+C3(i,i+1));
end
C1 = diag(C1);
    
C = C1 + C4 ;


l = 2.05;               % piston rod
r = 2.05/2;             % crankshaft, stroke/2
M_connecting_rod = 1390;  % (kg)
M_balance_weight_on_crankshaft = (460+575)/2;             % (kg), apo 460 ews 575
Mrot = M_connecting_rod + M_balance_weight_on_crankshaft;   % peristrefomeni maza
M_crosshead = 1464;                       % (kg)
M_piston_Mpiston_rod = 1011;               % (kg), piston+piston rod
Mrec = M_crosshead + M_piston_Mpiston_rod ;  % palindromiki maza 3330
B = 0.5;      % bore

% diafora fasis kilindron
n = 5;
aa = 360/n;
a = aa*pi/180;

Ds = Ds1(y(1:12),r,l);              % kalw to function Ds1
Ds_tonos = Ds1_tonos(y(1:12),r,l);  % kalw to function Ds1_tonos

VI = diag(Mrot*Ds.^2);    % Variable Inertia
CI = diag(0.5*Mrot*r^2);  % Constant Inertia
Q = VI + J;         % pinakas 16x16 , diagonios, epilego an thelo ypologismo me Variable h Constant inertia
%Q = J;  % diairw mono me J otan den epilego to Tinertia


% diadikasia gia tin entoli ode
main_func = zeros(24,1); % arxikopoiisi
main_func(1:12) = y(13:24);
%rpm_ = 60*mean(y(13:24))/(2*pi);
rpm_ = mean(y(14:21)); 
%rpm_ = y(13);
%% vasiki eksisosi limeni ws pros thita''

if time < 3
    
    main_func(13:24) = inv(Q) * ( - C*y(13:24) - K*y(1:12) - Tprop(y(24)) + Tgas25(y(1:12),Ds) - Tinertia(Mrec,Ds,Ds_tonos,y(13:24)) - Tfriction(y(13:24)) ) ;
    
else
    
    main_func(13:24) = inv(Q) * ( - C*y(13:24) - K*y(1:12) - Tprop(y(24)) + Tgas(y(1:12),Ds,rpm_) - Tinertia(Mrec,Ds,Ds_tonos,y(13:24)) - Tfriction(y(13:24)) ) ;
    
end

%if time<10
%    main_func(13:24) = inv(Q) * ( - C*y(13:24) - K*y(1:12) - Tprop(y(24)) + Tgas25(y(1:12),Ds) - Tinertia(Mrec,Ds,Ds_tonos,y(13:24)) - Tfriction(y(13:24)) ) ;
%else
%    main_func(13:24) = inv(Q) * ( - C*y(13:24) - K*y(1:12) - Tprop(y(24)) + Tgas75(y(1:12),Ds) - Tinertia(Mrec,Ds,Ds_tonos,y(13:24)) - Tfriction(y(13:24)) ) ;
%end  
%main_func(13:24) = inv(Q) * ( - C*y(13:24) - K*y(1:12) - Tprop(y(24)) + Tstart(y(1:12),Ds) + Tgas(y(1:12),Ds) - Tinertia(Mrec,Ds,Ds_tonos,y(13:24)) - Tfriction(y(13:24)) ) ;
%main_func(13:24) = inv(Q) * ( - C*y(13:24) - K*y(1:12) - Tprop(y(24)) + Tstart(y(1:12),Ds) - Tinertia(Mrec,Ds,Ds_tonos,y(13:24)) - Tfriction(y(13:24)) ) ;  % starting engine
% pinakas y  :  1:N gonies kathe DoF, N+1:2*N goniakes taxitites kathe DoF
end

