clear all
close all
clc

% N plithos DoF, N=12

% idioi pinakes inertia, stiffness, damping kai sto main_func (dedomena gia 5RTflex50D)
% pinakas y  :  sthles 1:N gonies kathe DoF, sthles N+1:2*N goniakes taxitites kathe DoF
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


% ypologismos idiosixnotinon systimatos
[eigenVectors, eigenValues] = eig(K,J);                       % eigenValues = lamda^2 (det(K-lamda^2*J)=0)
idiosixnotites = sqrt(eigenValues);                           % rad/s
idiosixnotites_rpm = 60*idiosixnotites/(2*pi);                % rpm
orizousa = det(K*eigenVectors - J*eigenVectors*eigenValues);  % epalithefsi oti vgainei miden h orizousa
idiosixnotites_rpm2 = diag(idiosixnotites_rpm);               % metatrepei ton diagonio pinaka se pinaka sthlh

dtmax = 10^(-4);    
opts = odeset('MaxStep',dtmax,'RelTol',1e-8,'AbsTol',1e-8);


initial_conditions = zeros(1,24)'; % arxikopoiisi

% arxikes sinthikes gonias (ta DoF idia arxiki sinthiki gonias)

initial_conditions(1:12) = 2*2*pi/360;               % arxiki sinthiki gonias (ola ta DoF idia arxiki sinthiki gonias)
%initial_conditions(1:12) = 37*2*pi/360;    % starting engine
% 37 moires gonia wste o kilindros 5 na vrisketai sto Ano Nekro Simeio, 
% gonia 37 moiron epilextike me dokimes me vasi to Ds tou 5ou kilindrou, (tyxaia epilogi o 5os kilindros)

start_rpm = 60;
initial_conditions(13:24) = 2*pi*start_rpm/60;    % arxiki sinthiki goniakis taxititas (ola ta DoF idia arxiki sinthiki gonias)
%initial_conditions(13:24) = 0;                   % starting engine

[time,y] = ode23tb(@main_func,[0 20],initial_conditions,opts); % epilisi vasikis eksisosis me ode
for i=1:length(y)
    y_main = main_func(time(i),y(i,:)');
    y_dotdot(i,:) = y_main(13:24);
end
%[time,y] = ode23tb(@main_func,[0 1.59],initial_conditions,opts);  %starting system, (time=1.65 sec)
% [0 20], 20 seconds : xronos pou trexei

figure
plot(time,y_dotdot(:,4))
xlabel('Time (s)')
ylabel('Anglular Accel (rad/s^2)')
title('Angular Accel - Time 4rth DOF')

figure
hold on
for i=1:size(y_dotdot,2)
    plot(time,y_dotdot(:,i));
end
hold off
xlabel('Time (s)')
ylabel('Anglular Accel (rad/s^2)')
title('Angular Accel - Time')
% plotarei ta stoixeia toy pinaka y(1:N), dhladh tis gonies se rad me to xrono
figure
hold on
% plot(time,y(:,1))
% plot(time,y(:,2))
% plot(time,y(:,3))
% plot(time,y(:,4))
% plot(time,y(:,5))
% plot(time,y(:,6))
% plot(time,y(:,7))
% plot(time,y(:,8))
% plot(time,y(:,9))
% plot(time,y(:,10))
% plot(time,y(:,11))
% plot(time,y(:,12))
% subplot(1,2,1)
% plot(time,y(:,1),time,y(:,2),time,y(:,3),time,y(:,4),time,y(:,5),time,y(:,6),time,y(:,7),time,y(:,8),time,y(:,9),time,y(:,10),time,y(:,11),time,y(:,12))
% 
% xlabel('Time (s)')
% ylabel('Angles (rad)')
% title('Angles - Time')

% plotarei ta stoixeia toy pinaka y(1:N), dhladh tis gonies se degrees me to xrono
% figure
% hold on
% plot(time,y(:,1)*360/(2*pi))
% plot(time,y(:,2)*360/(2*pi))
% plot(time,y(:,3)*360/(2*pi))
% plot(time,y(:,4)*360/(2*pi))
% plot(time,y(:,5)*360/(2*pi))
% plot(time,y(:,6)*360/(2*pi))
% plot(time,y(:,7)*360/(2*pi))
% plot(time,y(:,8)*360/(2*pi))
% plot(time,y(:,9)*360/(2*pi))
% plot(time,y(:,10)*360/(2*pi))
% plot(time,y(:,11)*360/(2*pi))
% plot(time,y(:,12)*360/(2*pi))

% subplot(1,2,2)
plot(time,y(:,1)*360/(2*pi),time,y(:,2)*360/(2*pi),time,y(:,3)*360/(2*pi),time,y(:,4)*360/(2*pi),time,y(:,5)*360/(2*pi),time,y(:,6)*360/(2*pi),time,y(:,7)*360/(2*pi),time,y(:,8)*360/(2*pi),time,y(:,9)*360/(2*pi),time,y(:,10)*360/(2*pi),time,y(:,11)*360/(2*pi),time,y(:,12)*360/(2*pi))

xlabel('Time (s)')
ylabel('Angles (degrees)')
title('Angles - Time')


% plotarei ta stoixeia toy pinaka y(N+1:2*N), dhladh tis goniakes taxitites se rad/s me to xrono
figure
% hold on
% plot(time,y(:,13))
% plot(time,y(:,14))
% plot(time,y(:,15))
% plot(time,y(:,16))
% plot(time,y(:,17))
% plot(time,y(:,18))
% plot(time,y(:,19))
% plot(time,y(:,20))
% plot(time,y(:,21))
% plot(time,y(:,22))
% plot(time,y(:,23))
% plot(time,y(:,24))

% subplot(2,2,1)
% plot(time,y(:,13),time,y(:,14),time,y(:,15),time,y(:,16),time,y(:,17),time,y(:,18),time,y(:,19),time,y(:,20),time,y(:,21),time,y(:,22),time,y(:,23),time,y(:,24))
% 
% xlabel('Time (s)')
% ylabel('Anglular Velocities (rad/s)')
% title('Angular Velocities - Time')

% SMOOTHDATA plotarei ta stoixeia toy pinaka y(N+1:2*N), dhladh tis goniakes taxitites se rad/s me to xrono
% figure
% hold on
y13_smooth = smoothdata(y(:,13));
y14_smooth = smoothdata(y(:,14));
y15_smooth = smoothdata(y(:,15));
y16_smooth = smoothdata(y(:,16));
y17_smooth = smoothdata(y(:,17));
y18_smooth = smoothdata(y(:,18));
y19_smooth = smoothdata(y(:,19));
y20_smooth = smoothdata(y(:,20));
y21_smooth = smoothdata(y(:,21));
y22_smooth = smoothdata(y(:,22));
y23_smooth = smoothdata(y(:,23));
y24_smooth = smoothdata(y(:,24));

% plot(time,y13_smooth)
% plot(time,y14_smooth)
% plot(time,y15_smooth)
% plot(time,y16_smooth)
% plot(time,y17_smooth)
% plot(time,y18_smooth)
% plot(time,y19_smooth)
% plot(time,y20_smooth)
% plot(time,y21_smooth)
% plot(time,y22_smooth)
% plot(time,y23_smooth)
% plot(time,y24_smooth)
% subplot(2,2,2)
% plot(time,y13_smooth,time,y14_smooth,time,y15_smooth,time,y16_smooth,time,y17_smooth,time,y18_smooth,time,y19_smooth,time,y20_smooth,time,y21_smooth,time,y22_smooth,time,y23_smooth,time,y24_smooth)
% 
% xlabel('Time (s)')
% ylabel('Anglular Velocities (rad/s)')
% title('Smooth : Angular Velocities - Time')

% plotarei ta stoixeia toy pinaka y(N+1:2*N), dhladh tis goniakes taxitites se rpm me to xrono
% figure
% hold on
% plot(time,y(:,13)*60/(2*pi))
% plot(time,y(:,14)*60/(2*pi))
% plot(time,y(:,15)*60/(2*pi))
% plot(time,y(:,16)*60/(2*pi))
% plot(time,y(:,17)*60/(2*pi))
% plot(time,y(:,18)*60/(2*pi))
% plot(time,y(:,19)*60/(2*pi))
% plot(time,y(:,20)*60/(2*pi))
% plot(time,y(:,21)*60/(2*pi))
% plot(time,y(:,22)*60/(2*pi))
% plot(time,y(:,23)*60/(2*pi))
% plot(time,y(:,24)*60/(2*pi))

%subplot(2,2,3)
plot(time,y(:,13)*60/(2*pi),time,y(:,14)*60/(2*pi),time,y(:,15)*60/(2*pi),time,y(:,16)*60/(2*pi),time,y(:,17)*60/(2*pi),time,y(:,18)*60/(2*pi),time,y(:,19)*60/(2*pi),time,y(:,20)*60/(2*pi),time,y(:,21)*60/(2*pi),time,y(:,22)*60/(2*pi),time,y(:,23)*60/(2*pi),time,y(:,24)*60/(2*pi))

xlabel('Time (s)')
ylabel('Anglular Velocities (rpm)')
title('Angular Velocities - Time')

% SMOOTHDATA plotarei ta stoixeia toy pinaka y(N+1:2*N), dhladh tis goniakes taxitites se rpm me to xrono
% figure
% hold on
% plot(time,y13_smooth*60/(2*pi))
% plot(time,y14_smooth*60/(2*pi))
% plot(time,y15_smooth*60/(2*pi))
% plot(time,y16_smooth*60/(2*pi))
% plot(time,y17_smooth*60/(2*pi))
% plot(time,y18_smooth*60/(2*pi))
% plot(time,y19_smooth*60/(2*pi))
% plot(time,y20_smooth*60/(2*pi))
% plot(time,y21_smooth*60/(2*pi))
% plot(time,y22_smooth*60/(2*pi))
% plot(time,y23_smooth*60/(2*pi))
% plot(time,y24_smooth*60/(2*pi))

% subplot(2,2,4)
% plot(time,y13_smooth*60/(2*pi),time,y14_smooth*60/(2*pi),time,y15_smooth*60/(2*pi),time,y16_smooth*60/(2*pi),time,y17_smooth*60/(2*pi),time,y18_smooth*60/(2*pi),time,y19_smooth*60/(2*pi),time,y20_smooth*60/(2*pi),time,y21_smooth*60/(2*pi),time,y22_smooth*60/(2*pi),time,y23_smooth*60/(2*pi),time,y24_smooth*60/(2*pi))
% 
% xlabel('Time (s)')
% ylabel('Anglular Velocities (rpm)')
% title('Smooth : Angular Velocities - Time')


% plotarei to Tgas kai tin Piesi me tin gonia gia opoio kilindro diallekso (oloi oi kilindroi)
figure
hold on
for k = 1:length(y)
    DsDs = Ds1(y(k,:),r,l);
    rpm_ = mean(y(k,14:24));
    A = Tgas(y(k,1:12),DsDs,rpm_);
    
    AA1(k) = A(3); % epilego poio Ds thelo gia plot
    PP1(k) = -AA1(k)/(DsDs(3,3)*pi/4*B^2); % piesi (gia plot)
    DSDS1(k) = -DsDs(3,3);
    
    AA2(k) = A(4); % epilego poio Ds thelo gia plot
    PP2(k) = -AA2(k)/(DsDs(4,4)*pi/4*B^2); % piesi (gia plot)
    DSDS2(k) = -DsDs(4,4);

    AA3(k) = A(5); % epilego poio Ds thelo gia plot
    PP3(k) = -AA3(k)/(DsDs(5,5)*pi/4*B^2); % piesi (gia plot)
    DSDS3(k) = -DsDs(5,5);
    
    AA4(k) = A(6); % epilego poio Ds thelo gia plot
    PP4(k) = -AA4(k)/(DsDs(6,6)*pi/4*B^2); % piesi (gia plot)
    DSDS4(k) = -DsDs(6,6);

    AA5(k) = A(7); % epilego poio Ds thelo gia plot
    PP5(k) = -AA5(k)/(DsDs(7,7)*pi/4*B^2); % piesi (gia plot)
    DSDS5(k) = -DsDs(7,7);

end

% plot(y(:,3),AA1)
% plot(y(:,4),AA2)
% plot(y(:,5),AA3)
% plot(y(:,6),AA4)
% plot(y(:,7),AA5)

subplot(1,2,1)
plot(y(:,3),AA1,y(:,4),AA2,y(:,5),AA3,y(:,6),AA4,y(:,7),AA5)

xlabel('Angle (rad)')
ylabel('Gas Torques (Nm)')
title('Gas Torques - Angle')

% figure
% hold on
% plot(y(:,3),PP1)
% plot(y(:,4),PP2)
% plot(y(:,5),PP3)
% plot(y(:,6),PP4)
% plot(y(:,7),PP5)

subplot(1,2,2)
plot(y(:,3),PP1,y(:,4),PP2,y(:,5),PP3,y(:,6),PP4,y(:,7),PP5)

xlabel('Angle (rad)')
ylabel('Cylinder Pressures (Pa)')
title('Cylinder Pressures - Angle')

% plotarei to Tgas kai tin Piesi me ton xrono gia opoio kilindro diallekso (oloi oi kilindroi)
% figure
% hold on
% plot(time,AA1)
% plot(time,AA2)
% plot(time,AA3)
% plot(time,AA4)
% plot(time,AA5)

% subplot(2,2,3)
% plot(time,AA1,time,AA2,time,AA3,time,AA4,time,AA5)
% 
% xlabel('Time (s)')
% ylabel('Gas Torques (Nm)')
% title('Gas Torques - Time')

% figure
% hold on
% plot(time,PP1)
% plot(time,PP2)
% plot(time,PP3)
% plot(time,PP4)
% plot(time,PP5)

% subplot(2,2,4)
% plot(time,PP1,time,PP2,time,PP3,time,PP4,time,PP5)
% 
% xlabel('Time (s)')
% ylabel('Cylinder Pressures (Pa)')
% title('Cylinder Pressures - Time')


% plotarei tin Tprop (ropi elikas) me tin goniaki taxitita
figure
for p = 1:length(y)
    Tp = Tprop(y(p,24));
    TP(p) = Tp(12);   % teleutaio DoF pou se auto askeitai to Tprop
end

% subplot(1,2,1)
% plot(y(:,24),TP)
% xlabel('Angular Velocity (rad/s)')
% ylabel('Propeller Torque (Nm)')
% title('Propeller Torque - Angular Velocity')

%figure
%subplot(1,2,2)
plot(60*y(:,24)/(2*pi),TP)
xlabel('Angular Velocity (rpm)')
ylabel('Propeller Torque (Nm)')
title('Propeller Torque - Angular Velocity')

% xreiazontai pio kato
PMCR = 6350*10^3;
NMCR = 99;
WMCR = 2*pi*NMCR/60;
vima = 0.1 ;

% plotarei tin Tprop me ton xrono + SMOOTHDATA plot
figure
hold on
plot(time,TP)
Tprop_smooth = smoothdata(TP);

plot(time,Tprop_smooth)
xlabel('Time (s)')
ylabel('Propeller Torque')
legend({'Propeller Torque','Smooth : Propeller Torque'},'Location','southeast')
title('Propeller Torque - Time')


% plotarei to Tfriction gia opoio stoixeio diallekso sinartisi tis goniakis taxititas se rad/s kai se rpm
figure
hold on
for k = 1:length(y)
    
    F1 = Tfriction(y(k,13));
    Fr1(k) = F1(1,1);  % epilego se poio DoF plotaro to Tfriction
    
    F2 = Tfriction(y(k,14));
    Fr2(k) = F2(2,2);
    
    F3 = Tfriction(y(k,15));
    Fr3(k) = F3(3,3);
    
    F4 = Tfriction(y(k,16));
    Fr4(k) = F4(4,4);
    
    F5 = Tfriction(y(k,17));
    Fr5(k) = F5(5,5);
    
    F6 = Tfriction(y(k,18));
    Fr6(k) = F6(6,6);
    
    F7 = Tfriction(y(k,19));
    Fr7(k) = F7(7,7);
    
    F8 = Tfriction(y(k,20));
    Fr8(k) = F8(8,8);
    
    F9 = Tfriction(y(k,21));
    Fr9(k) = F9(9,9);
    
    F10 = Tfriction(y(k,22));
    Fr10(k) = F10(10,10);
    
    F11 = Tfriction(y(k,23));
    Fr11(k) = F11(11,11);
    
    F12 = Tfriction(y(k,24));
    Fr12(k) = F12(12,12);
end
% plot(y(:,13),Fr1)
% plot(y(:,14),Fr2)
% plot(y(:,15),Fr3)
% plot(y(:,16),Fr4)
% plot(y(:,17),Fr5)
% plot(y(:,18),Fr6)
% plot(y(:,19),Fr7)
% plot(y(:,20),Fr8)
% plot(y(:,21),Fr9)
% plot(y(:,22),Fr10)
% plot(y(:,23),Fr11)
% plot(y(:,24),Fr12)

subplot(1,2,1)
plot(y(:,13),Fr1,y(:,14),Fr2,y(:,15),Fr3,y(:,16),Fr4,y(:,17),Fr5,y(:,18),Fr6,y(:,19),Fr7,y(:,20),Fr8,y(:,21),Fr9,y(:,22),Fr10,y(:,23),Fr11,y(:,24),Fr12)

xlabel('Angular Velocity (rad/s)')
ylabel('Friction Torques (Nm)')
title ('Friction Torques - Angular Velocity')

% figure
% hold on
% plot(y(:,13)*60/(2*pi),Fr1)
% plot(y(:,14)*60/(2*pi),Fr2)
% plot(y(:,15)*60/(2*pi),Fr3)
% plot(y(:,16)*60/(2*pi),Fr4)
% plot(y(:,17)*60/(2*pi),Fr5)
% plot(y(:,18)*60/(2*pi),Fr6)
% plot(y(:,19)*60/(2*pi),Fr7)
% plot(y(:,20)*60/(2*pi),Fr8)
% plot(y(:,21)*60/(2*pi),Fr9)
% plot(y(:,22)*60/(2*pi),Fr10)
% plot(y(:,23)*60/(2*pi),Fr11)
% plot(y(:,24)*60/(2*pi),Fr12)

subplot(1,2,2)
plot(y(:,13)*60/(2*pi),Fr1,y(:,14)*60/(2*pi),Fr2,y(:,15)*60/(2*pi),Fr3,y(:,16)*60/(2*pi),Fr4,y(:,17)*60/(2*pi),Fr5,y(:,18)*60/(2*pi),Fr6,y(:,19)*60/(2*pi),Fr7,y(:,20)*60/(2*pi),Fr8,y(:,21)*60/(2*pi),Fr9,y(:,22)*60/(2*pi),Fr10,y(:,23)*60/(2*pi),Fr11,y(:,24)*60/(2*pi),Fr12)

xlabel('Angular Velocity (rpm)')
ylabel('Friction Torques (Nm)')
title ('Friction Torques - Angular Velocity')

% plotarei to Tfriction gia opoio stoixeio diallekso sinartisi tou xronou
figure
hold on
plot(time,Fr1)
plot(time,Fr2)
plot(time,Fr3)
plot(time,Fr4)
plot(time,Fr5)
plot(time,Fr6)
plot(time,Fr7)
plot(time,Fr8)
plot(time,Fr9)
plot(time,Fr10)
plot(time,Fr11)
plot(time,Fr12)

% SMOOTHDATA plotarei to Tfriction gia opoio stoixeio diallekso sinartisi tou xronou
fr13 = smoothdata(Fr1);
plot(time,fr13)
fr14 = smoothdata(Fr2);
plot(time,fr14)
fr15 = smoothdata(Fr3);
plot(time,fr15)
fr16 = smoothdata(Fr4);
plot(time,fr16)
fr17 = smoothdata(Fr5);
plot(time,fr17)
fr18 = smoothdata(Fr6);
plot(time,fr18)
fr19 = smoothdata(Fr7);
plot(time,fr19)
fr20 = smoothdata(Fr8);
plot(time,fr20)
fr21 = smoothdata(Fr9);
plot(time,fr21)
fr22 = smoothdata(Fr10);
plot(time,fr22)
fr23 = smoothdata(Fr11);
plot(time,fr23)
fr24 = smoothdata(Fr12);
plot(time,fr24)

xlabel('Time (s)')
ylabel('Friction Torque (Nm)')
%legend({'Smooth : Friction Torque','Friction Torque'},'Location','southeast')
title ('Friction Torque - Time')

figure
hold on
plot(time,fr15,'k')
plot(time,Fr3,'r')
xlabel('Time (s)')
ylabel('Friction Torque (Nm)')
legend({'Smooth : Friction Torque','Friction Torque'},'Location','southeast')
title ('Friction Torque - Time')


% plotarei to Tinertia gia opoio stoixeio diallekso sinartisi tou xronou
figure
hold on
for n = 1:length(y)
    DsDs = Ds1(y(n,:),r,l);
    DsDstonos = Ds1_tonos(y(n,:),r,l);
    In = Tinertia(Mrec,DsDs,DsDstonos,y(n,13:24)');
    
    Inertia1(n) = In(3);
    Inertia2(n) = In(4);
    Inertia3(n) = In(5);
    Inertia4(n) = In(6);
    Inertia5(n) = In(7);
end
subplot(1,2,1)
plot(time,Inertia1,time,Inertia2,time,Inertia3,time,Inertia4,time,Inertia5)
% plot(time,Inertia1)
% plot(time,Inertia2)
% plot(time,Inertia3)
% plot(time,Inertia4)
% plot(time,Inertia5)

xlabel('Time (s)')
ylabel('Reciprocating Masses Inertia Torques (Nm)')
title ('Reciprocating Masses Inertia Torques - Time')

%figure
subplot(1,2,2)
plot(time,Inertia1)
xlabel('Time (s)')
ylabel('Reciprocating Mass Inertia Torque (Nm)')
title ('Reciprocating Mass Inertia Torque (1 cylinder) - Time')

% piesi kai ropi ekkinisis logo pepiesmenou aera
figure
hold on
for kk = 1:length(y)
    DsDs = Ds1(y(kk,:),r,l);
    Astart = Tstart(y(kk,1:12),DsDs);
    
    AAstart(kk) = Astart(5); % epilego poio Ds thelo gia plot
    PPstart(kk) = -AAstart(kk)/(DsDs(5,5)*pi/4*B^2); % piesi (gia plot)

end

subplot(1,2,1)
plot(y(:,5),AAstart)
xlabel('Angle (rad)')
ylabel('Air Starting Torque (Nm)')
title('Air Starting Torque - Angle')

%figure
subplot(1,2,2)
plot(y(:,5),PPstart)
xlabel('Angle (rad)')
ylabel('Air Starting Pressure (Pa)')
title('Air Starting Pressure - Angle')

%figure
% subplot(2,2,3)
% plot(time,AAstart)
% 
% xlabel('Time (s)')
% ylabel('Air Starting Torque (Nm)')
% title('Air Starting Torque - Time')

%figure
% subplot(2,2,4)
% plot(time,PPstart)
% 
% xlabel('Time (s)')
% ylabel('Air Starting Pressure (Pa)')
% title('Air Starting Pressure - Time')
% 
% title('Air Starting Pressure - Time')



% plotarei tin piesi enos kilindrou sinartisi tis gonias
figure
subplot(1,2,1)
plot(y(:,3),PP1)

xlabel('Angle (rad)')
ylabel('Cylinder Pressures (Pa)')
title('Cylinder Pressure (1 cylinder) - Angle')


% plotarei to Tgas enos kilindrou sinartisi tis gonias
%figure
subplot(1,2,2)
plot(y(:,3),AA1)

xlabel('Angle (rad)')
ylabel('Gas Torque (Nm)')
title('Gas Torques (1 cylinder) - Angle')


% plotarei tin piesi enos kilindrou sinartisi tou xronou 
%figure
% subplot(2,2,3)
% plot(time,PP1)
% xlabel('Time (s)')
% ylabel('Cylinder Pressure (Pa)')
% title('Cylinder Pressure (1 cylinder) - Time')


% plotarei to Tgas enos kilindrou sinartisi tou xronou 
%figure
% subplot(2,2,4)
% plot(time,AA1)
% %title('Tgas for 1 cylinder - time')
% 
% xlabel('Time (s)')
% ylabel('Gas Torque (Nm)')
% title('Gas Torque (1 cylinder) - Time')

% Loading Diagram
C11 = 0.5 ;
C22 = 0.5 ;
C21 = 0.3 ;
C32 = 1.111 ;
C31 = -0.067 ;
C41 = 1 ;
C50 = 1 ;

vima = 0.5 ;

n1 = [0:vima:0.4*NMCR];
P1 = PMCR*C11*n1/NMCR;

n2 = [0.4*NMCR:vima:0.6*NMCR];
P2 = PMCR*(C22*n2.^2/NMCR^2+C21*n2/NMCR);

n3 = [0.6*NMCR:vima:0.96*NMCR];
P3 = PMCR*(C32*n3.^2/NMCR^2+C31*n3/NMCR);

n4 = [0.96*NMCR:vima:1*NMCR];
P4 = PMCR*C41*n4/NMCR;

n5 = [1*NMCR:vima:1.04*NMCR];
P5 = [C50*PMCR]*ones(size(n5));

P6 = [0:2000:1*PMCR];
n6 = 1.04*NMCR*ones(size(P6));


% overload
C11over = 0.6 ;
C22over = 0.33 ;
C21over = 0.468 ;
C32over = 1.111 ;
C41over = 1.066 ;
C50over = 1.1 ;

n1over = [0:0.1:0.4*NMCR];
P1over = PMCR*C11over*n1over/NMCR;

n2over = [0.4*NMCR:0.1:0.6*NMCR];
P2over = PMCR*(C22over*n2over.^2/NMCR^2+C21over*n2over/NMCR);

n3over = [0.6*NMCR:0.1:0.96*NMCR];
P3over = PMCR*(C32over*n3over.^2/NMCR^2);

n4over = [0.96*NMCR:0.1:1.032*NMCR];
P4over = PMCR*C41over*n4over/NMCR;

n5over = [1.032*NMCR:0.1:1.08*NMCR];
P5over = [C50over*PMCR]*ones(size(n5over));

P6over = [0:2000:1.1*PMCR];
n6over = 1.08*NMCR*ones(size(P6over));

figure
hold on

% plot Loading Diagram
ld(1) = plot(n1,P1,'b')
ld(2) = plot(n2,P2,'b')
ld(3) = plot(n3,P3,'b')
ld(4) = plot(n4,P4,'b')
ld(5) = plot(n5,P5,'b')
ld(6) = plot(n6,P6,'b')

% plot Loading Diagram, overload
ld(7) = plot(n1over,P1over,'k')
ld(8) = plot(n2over,P2over,'k')
ld(9) = plot(n3over,P3over,'k')
ld(10) = plot(n4over,P4over,'k')
ld(11) = plot(n5over,P5over,'k')
ld(12) = plot(n6over,P6over,'k')

% plotarei smooth tin ropi tis elikas
Power_propeller = y24_smooth.*Tprop_smooth';
ld(13) = plot(60*y24_smooth/(2*pi),Power_propeller,'r');

% plotarei smooth tin ropi tou kinitira, Pengine = Pgas - Pfriction
Tengine_start = AAstart;
Tengine_start_smooth = smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(Tengine_start)))))))))))))))));
Tengine = AA1+AA2+AA3+AA4+AA5;

Tengine_smooth_2 = smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(Tengine))))))))))))))))))))))))))))))));
%Tengine_smooth = smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(smoothdata(Tengine_smooth_2)))))))))))))))))))))))))))))));
Tengine_smooth =smoothdata(Tengine); 
Tfriction = Fr1+Fr2+Fr3+Fr4+Fr5+Fr6+Fr7+Fr8+Fr9+Fr10+Fr11+Fr12;
Tfriction_smooth = smoothdata(Tfriction);

ww = y21_smooth; % y21 flywheel (teliko simeio tis mixanis, meta ksekinaei to aksoniko
%ww = smoothdata(y21_smooth); % gia metavoli Load apo 75% se 100%

Power_engine = (Tengine_smooth - Tfriction_smooth).*ww';

Power_engine_start = (Tengine_start_smooth - Tfriction_smooth).*y21_smooth';  % starting system


ld(14) = plot(60*ww/(2*pi),Power_engine,'g');
%plot(60*y20_smooth/(2*pi),Power_engine_start,'c')

xlabel('Angular Velocity (rpm)')
ylabel('Power (Watt)')
legend(ld([1 7 13 14]),'Loading Diagram','Overload Diagram','Propeller Law','Engine Power','Location','northwest')
title('Loading Diagram - Propeller Law - Engine Power')

% figure
% hold on
% plot(time,Tengine_smooth,'b')
% xlabel('Time (s)')
% ylabel('Engine Torque (Nm)')
% title('Smooth Engine Torque - Time')

RPM = y24_smooth*60/(2*pi);

% plotarei to Ds
% figure
% hold on
% % plot(y(:,3),DSDS1,'k')
% % plot(y(:,4),DSDS2,'b')
% % plot(y(:,5),DSDS3,'r')
% % plot(y(:,6),DSDS4,'g')
% % plot(y(:,7),DSDS5,'c')
% 
% subplot(1,2,1)
% plot(y(:,3),DSDS1,'k',y(:,4),DSDS2,'b',y(:,5),DSDS3,'r',y(:,6),DSDS4,'g',y(:,7),DSDS5,'c')
% 
% xlabel('Angle (rad)')
% ylabel('Ds')
% title('Ds - Angle')
% 
% %figure
% subplot(1,2,2)
% plot(y(:,3),DSDS1)
% xlabel('Angle (rad)')
% ylabel('Ds')
% title('Ds (1 cylinder) - Angle')


figure
subplot(3,4,1)
% plot(time,y(:,13)*60/(2*pi),time,y13_smooth*60/(2*pi))
plot(time,y(:,13)*60/(2*pi))
% plot(time,y13_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Damper')

%figure
subplot(3,4,2)
% plot(time,y(:,14)*60/(2*pi),time,y14_smooth*60/(2*pi))
plot(time,y(:,14)*60/(2*pi))
% plot(time,y14_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Flange of M/E')

%figure
subplot(3,4,3)
% plot(time,y(:,15)*60/(2*pi),time,y15_smooth*60/(2*pi))
plot(time,y(:,15)*60/(2*pi))
% plot(time,y15_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Cylinder 1')

%figure
subplot(3,4,4)
% plot(time,y(:,16)*60/(2*pi),time,y16_smooth*60/(2*pi))
plot(time,y(:,16)*60/(2*pi))
% plot(time,y16_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Cylinder 2')

%figure
subplot(3,4,5)
% plot(time,y(:,17)*60/(2*pi),time,y17_smooth*60/(2*pi))
plot(time,y(:,17)*60/(2*pi))
% plot(time,y17_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Cylinder 3')

%figure
subplot(3,4,6)
% plot(time,y(:,18)*60/(2*pi),time,y18_smooth*60/(2*pi))
plot(time,y(:,18)*60/(2*pi))
% plot(time,y18_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Cylinder 4')

%figure
subplot(3,4,7)
% plot(time,y(:,19)*60/(2*pi),time,y19_smooth*60/(2*pi))
plot(time,y(:,19)*60/(2*pi))
% plot(time,y19_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Cylinder 5')

%figure
subplot(3,4,8)
% plot(time,y(:,20)*60/(2*pi),time,y20_smooth*60/(2*pi))
plot(time,y(:,20)*60/(2*pi))
% plot(time,y20_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Thrust Bearing')

%figure
subplot(3,4,9)
% plot(time,y(:,21)*60/(2*pi),time,y21_smooth*60/(2*pi))
plot(time,y(:,21)*60/(2*pi))
% plot(time,y21_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Flywheel')

%figure
subplot(3,4,10)
% plot(time,y(:,22)*60/(2*pi),time,y22_smooth*60/(2*pi))
plot(time,y(:,22)*60/(2*pi))
% plot(time,y22_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Intermidiate Shaft')

%figure
subplot(3,4,11)
% plot(time,y(:,23)*60/(2*pi),time,y23_smooth*60/(2*pi))
plot(time,y(:,23)*60/(2*pi))
% plot(time,y23_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('Propeller Shaft')

%figure
subplot(3,4,12)
% plot(time,y(:,24)*60/(2*pi),time,y24_smooth*60/(2*pi))
plot(time,y(:,24)*60/(2*pi))
% plot(time,y24_smooth*60/(2*pi))
xlabel('Time')
ylabel('Angular Velocity (rpm)')
title('FP Propeller')


figure
hold on
subplot(1,2,1)
U1 = y(:,15).*DSDS1';
U2 = y(:,16).*DSDS2';
U3 = y(:,17).*DSDS3';
U4 = y(:,18).*DSDS4';
U5 = y(:,19).*DSDS5';
% plot(time,U1)
% plot(time,U2)
% plot(time,U3)
% plot(time,U4)
% plot(time,U5)
plot(time,U1,time,U2,time,U3,time,U4,time,U5)
xlabel('Time (s)')
ylabel('Piston Speed (m/s)')
title('Piston Speed - Time')

% figure
% hold on
% plot(y(:,3),U1)
% plot(y(:,4),U2)
% plot(y(:,5),U3)
% plot(y(:,6),U4)
% plot(y(:,7),U5)
% subplot(2,2,2)
% plot(y(:,3),U1,y(:,4),U2,y(:,5),U3,y(:,6),U4,y(:,7),U5)
% xlabel('Angle (degrees)')
% ylabel('Piston Speed (m/s)')
% title('Piston Speed - Angle')

% figure
subplot(1,2,2)
plot(time,U1)
xlabel('Time (s)')
ylabel('Piston Speed (m/s)')
title('Piston Speed (1 cylinder) - Time')

%figure
% subplot(2,2,4)
% plot(y(:,3),U1)
% xlabel('Angle (degrees)')
% ylabel('Piston Speed (m/s)')
% title('Piston Speed (1 cylinder) - Angle')


figure
hold on
%subplot(1,2,1)
C1 = (y13_smooth*60/(2*pi))*2*r/30;
C2 = (y14_smooth*60/(2*pi))*2*r/30;
C3 = (y15_smooth*60/(2*pi))*2*r/30;
C4 = (y16_smooth*60/(2*pi))*2*r/30;
C5 = (y17_smooth*60/(2*pi))*2*r/30;
% plot(time,C1)
% plot(time,C2)
% plot(time,C3)
% plot(time,C4)
% plot(time,C5)
plot(time,C1,time,C2,time,C3,time,C4,time,C5)
xlabel('Time (s)')
ylabel('Mean Piston Speed (m/s)')
title('Mean Piston Speed - Time')

% figure
% hold on
% plot(y(:,3),C1)
% plot(y(:,4),C2)
% plot(y(:,5),C3)
% plot(y(:,6),C4)
% plot(y(:,7),C5)
% subplot(1,2,2)
% plot(y(:,3),C1,y(:,4),C2,y(:,5),C3,y(:,6),C4,y(:,7),C5)
% xlabel('Angle (degrees)')
% ylabel('Mean Piston Speed (m/s)')
% title('Mean Piston Speed - Angle')

figure
hold on
tfly1 = K(9,10)*(y(:,10)-y(:,9));
tfly2 = K(8,9)*(y(:,9)-y(:,8));
plot(time,tfly1,'g')
plot(time,tfly2,'r')
xlabel('Time (s)')
ylabel('Flywheel Torque (Nm)')
legend({'Torque after Flywheel','Torque before Flywheel'},'Location','southeast')
title('Flywheel Torque - Time')


% connecting rod angle
b3 = asin(r/l*sin(y(:,3)+0*a));
b4 = asin(r/l*sin(y(:,4)+3*a));
b5 = asin(r/l*sin(y(:,5)+2*a));
b6 = asin(r/l*sin(y(:,6)+1*a));
b7 = asin(r/l*sin(y(:,7)+4*a));

% figure
% hold on
% subplot(1,2,1)
% % plot(time,b3*360/(2*pi))
% % plot(time,b4*360/(2*pi))
% % plot(time,b5*360/(2*pi))
% % plot(time,b6*360/(2*pi))
% % plot(time,b7*360/(2*pi))
% plot(time,b3*360/(2*pi),time,b4*360/(2*pi),time,b5*360/(2*pi),time,b6*360/(2*pi),time,b7*360/(2*pi))
% xlabel('Time (s)')
% ylabel('Connecting Rod Angles (degrees)')
% title('Connecting Rod Angles - Time')
% %figure
% subplot(1,2,2)
% plot(time,b3*360/(2*pi))
% xlabel('Time (s)')
% ylabel('Connecting Rod Angle (degrees)')
% title('Connecting Rod Angle (1 cylinder) - Time')


% gas force Fg
Fg1 = pi/4*B^2*PP1;
Fg2 = pi/4*B^2*PP2;
Fg3 = pi/4*B^2*PP3;
Fg4 = pi/4*B^2*PP4;
Fg5 = pi/4*B^2*PP5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Fg1)
% % plot(y(:,4),Fg2)
% % plot(y(:,5),Fg3)
% % plot(y(:,6),Fg4)
% % plot(y(:,7),Fg5)
% plot(y(:,3),Fg1,y(:,4),Fg2,y(:,5),Fg3,y(:,6),Fg4,y(:,7),Fg5)
% xlabel('Angle (degrees)')
% ylabel('Piston Gas Forces (N)')
% title('Piston Gas Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Fg1)
% xlabel('Angle (degrees)')
% ylabel('Piston Gas Force (N)')
% title('Piston Gas Force (1 cylinder) - Angle')


% Gas Lateral Guide Force FgN
FgN1 = Fg1.*tan(b3');
FgN2 = Fg2.*tan(b4');
FgN3 = Fg3.*tan(b5');
FgN4 = Fg4.*tan(b6');
FgN5 = Fg5.*tan(b7');

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),-FgN1)
% % plot(y(:,4),-FgN2)
% % plot(y(:,5),-FgN3)
% % plot(y(:,6),-FgN4)
% % plot(y(:,7),-FgN5)
% plot(y(:,3),-FgN1,y(:,4),-FgN2,y(:,5),-FgN3,y(:,6),-FgN4,y(:,7),-FgN5)
% xlabel('Angle (degrees)')
% ylabel('Gas Lateral Guide Forces (N)')
% title('Gas Lateral Guide Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),-FgN1)
% xlabel('Angle (degrees)')
% ylabel('Gas Lateral Guide Force (N)')
% title('Gas Lateral Guide Force (1 cylinder) - Angle')


% Gas Crankpin Force Fgs
Fgs1 = Fg1./cos(b3');
Fgs2 = Fg2./cos(b4');
Fgs3 = Fg3./cos(b5');
Fgs4 = Fg4./cos(b6');
Fgs5 = Fg5./cos(b7');

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Fgs1)
% % plot(y(:,4),Fgs2)
% % plot(y(:,5),Fgs3)
% % plot(y(:,6),Fgs4)
% % plot(y(:,7),Fgs5)
% plot(y(:,3),Fgs1,y(:,4),Fgs2,y(:,5),Fgs3,y(:,6),Fgs4,y(:,7),Fgs5)
% xlabel('Angle (degrees)')
% ylabel('Gas Crankpin Forces (N)')
% title('Gas Crankpin Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Fgs1)
% xlabel('Angle (degrees)')
% ylabel('Gas Crankpin Force (N)')
% title('Gas Crankpin Force (1 cylinder) - Angle')


% Gas Tangential Crankpin Force Fgt
Fgt1 = (1/r)*Fg1.*DSDS1;
Fgt2 = (1/r)*Fg2.*DSDS2;
Fgt3 = (1/r)*Fg3.*DSDS3;
Fgt4 = (1/r)*Fg4.*DSDS4;
Fgt5 = (1/r)*Fg5.*DSDS5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Fgt1)
% % plot(y(:,4),Fgt2)
% % plot(y(:,5),Fgt3)
% % plot(y(:,6),Fgt4)
% % plot(y(:,7),Fgt5)
% plot(y(:,3),Fgt1,y(:,4),Fgt2,y(:,5),Fgt3,y(:,6),Fgt4,y(:,7),Fgt5)
% xlabel('Angle (degrees)')
% ylabel('Gas Tangential Crankpin Forces (N)')
% title('Gas Tangential Crankpin Forces - Angle')
% %figure 
% subplot(1,2,2)
% plot(y(:,3),Fgt1)
% xlabel('Angle (degrees)')
% ylabel('Gas Tangential Crankpin Force (N)')
% title('Gas Tangential Crankpin Force (1 cylinder) - Angle')


% Gas Radial Crankpin Force FgR
FgR1 = -Fg1.*(cos(y(:,3)+0*a+b3)./cos(b3))';
FgR2 = -Fg2.*(cos(y(:,4)+3*a+b4)./cos(b4))';
FgR3 = -Fg3.*(cos(y(:,5)+2*a+b5)./cos(b5))';
FgR4 = -Fg4.*(cos(y(:,6)+1*a+b6)./cos(b6))';
FgR5 = -Fg5.*(cos(y(:,7)+4*a+b7)./cos(b7))';

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),FgR1)
% % plot(y(:,4),FgR2)
% % plot(y(:,5),FgR3)
% % plot(y(:,6),FgR4)
% % plot(y(:,7),FgR5)
% plot(y(:,3),FgR1,y(:,4),FgR2,y(:,5),FgR3,y(:,6),FgR4,y(:,7),FgR5)
% xlabel('Angle (degrees)')
% ylabel('Gas Radial Crankpin Forces (N)')
% title('Gas Radial Crankpin Forces - Angle')
% %figure 
% subplot(1,2,2)
% plot(y(:,3),FgR1)
% xlabel('Angle (degrees)')
% ylabel('Gas Radial Crankpin Force (N)')
% title('Gas Radial Crankpin Force (1 cylinder) - Angle')



% Tangential Reciprocating Mass Inertia Force Flt
Flt1 = Inertia1/r;
Flt2 = Inertia2/r;
Flt3 = Inertia3/r;
Flt4 = Inertia4/r;
Flt5 = Inertia5/r;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Flt1)
% % plot(y(:,4),Flt2)
% % plot(y(:,5),Flt3)
% % plot(y(:,6),Flt4)
% % plot(y(:,7),Flt5)
% plot(y(:,3),Flt1,y(:,4),Flt2,y(:,5),Flt3,y(:,6),Flt4,y(:,7),Flt5)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Masses Inertia Tangential Crankpin Forces (N)')
% title('Reciprocating Masses Inertia Tangential Crankpin Forces - Angle')
% %figure 
% subplot(1,2,2)
% plot(y(:,3),Flt1)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Mass Inertia Tangential Crankpin Force (N)')
% title('Reciprocating Mass Inertia Tangential Crankpin Force (1 cylinder) - Angle')


% Reciprocating Masses Inertia Forces Fl
Fl1 = r*Flt1./DSDS1;
Fl2 = r*Flt2./DSDS2;
Fl3 = r*Flt3./DSDS3;
Fl4 = r*Flt4./DSDS4;
Fl5 = r*Flt5./DSDS5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Fl1)
% % plot(y(:,4),Fl2)
% % plot(y(:,5),Fl3)
% % plot(y(:,6),Fl4)
% % plot(y(:,7),Fl5)
% plot(y(:,3),Fl1,y(:,4),Fl2,y(:,5),Fl3,y(:,6),Fl4,y(:,7),Fl5)
% xlabel('Angle (degrees)')
% ylabel('Piston Reciprocating Masses Inertia Forces (N)')
% title('Piston Reciprocating Masses Inertia Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Fl1)
% xlabel('Angle (degrees)')
% ylabel('Piston Reciprocating Mass Inertia Force (N)')
% title('Piston Reciprocating Mass Inertia Force (1 cylinder) - Angle')


% Reciprocating Masses Inertia Crankpin Forces Fls
Fls1 = Fl1./cos(b3');
Fls2 = Fl2./cos(b4');
Fls3 = Fl3./cos(b5');
Fls4 = Fl4./cos(b6');
Fls5 = Fl5./cos(b7');

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Fls1)
% % plot(y(:,4),Fls2)
% % plot(y(:,5),Fls3)
% % plot(y(:,6),Fls4)
% % plot(y(:,7),Fls5)
% plot(y(:,3),Fls1,y(:,4),Fls2,y(:,5),Fls3,y(:,6),Fls4,y(:,7),Fls5)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Masses Inertia Crankpin Forces (N)')
% title('Reciprocating Masses Inertia Crankpin Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Fls1)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Mass Inertia Crankpin Force (N)')
% title('Reciprocating Mass Inertia Crankpin Force (1 cylinder) - Angle')


% Radial Reciprocating Masses Inertia Forces FlR
FlR1 = -Fl1.*(cos(y(:,3)+0*a+b3)./cos(b3))';
FlR2 = -Fl2.*(cos(y(:,4)+3*a+b4)./cos(b4))';
FlR3 = -Fl3.*(cos(y(:,5)+2*a+b5)./cos(b5))';
FlR4 = -Fl4.*(cos(y(:,6)+1*a+b6)./cos(b6))';
FlR5 = -Fl5.*(cos(y(:,7)+4*a+b7)./cos(b7))';

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),FlR1)
% % plot(y(:,4),FlR2)
% % plot(y(:,5),FlR3)
% % plot(y(:,6),FlR4)
% % plot(y(:,7),FlR5)
% plot(y(:,3),FlR1,y(:,4),FlR2,y(:,5),FlR3,y(:,6),FlR4,y(:,7),FlR5)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Masses Inertia Radial Crankpin Forces (N)')
% title('Reciprocating Masses Inertia Radial Crankpin Forces - Angle')
% %figure 
% subplot(1,2,2)
% plot(y(:,3),FlR1)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Mass Inertia Radial Crankpin Force (N)')
% title('Reciprocating Mass Inertia Radial Crankpin Force (1 cylinder) - Angle')


%  Reciprocating Masses Inertia Lateral Guide Forces FlN
FlN1 = Fl1.*tan(b3');
FlN2 = Fl2.*tan(b4');
FlN3 = Fl3.*tan(b5');
FlN4 = Fl4.*tan(b6');
FlN5 = Fl5.*tan(b7');

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),-FlN1)
% % plot(y(:,4),-FlN2)
% % plot(y(:,5),-FlN3)
% % plot(y(:,6),-FlN4)
% % plot(y(:,7),-FlN5)
% plot(y(:,3),-FlN1,y(:,4),-FlN2,y(:,5),-FlN3,y(:,6),-FlN4,y(:,7),-FlN5)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Masses Inertia Lateral Guide Forces (N)')
% title('Reciprocating Masses Inertia Lateral Guide Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),-FlN1)
% xlabel('Angle (degrees)')
% ylabel('Reciprocating Mass Inertia Lateral Guide Force (N)')
% title('Reciprocating Mass Inertia Lateral Guide Force (1 cylinder) - Angle')


% total Force
Ftot1 = Fl1 + Fg1;
Ftot2 = Fl2 + Fg2;
Ftot3 = Fl3 + Fg3;
Ftot4 = Fl4 + Fg4;
Ftot5 = Fl5 + Fg5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Ftot1)
% % plot(y(:,4),Ftot2)
% % plot(y(:,5),Ftot3)
% % plot(y(:,6),Ftot4)
% % plot(y(:,7),Ftot5)
% plot(y(:,3),Ftot1,y(:,4),Ftot2,y(:,5),Ftot3,y(:,6),Ftot4,y(:,7),Ftot5)
% xlabel('Angle (degrees)')
% ylabel('Total Piston Forces (N)')
% title('Total Piston Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Ftot1)
% xlabel('Angle (degrees)')
% ylabel('Total Piston Force (N)')
% title('Total Piston Force (1 cylinder) - Angle')


% total lateral guide force
FtotN1 = (FlN1 + FgN1);
FtotN2 = (FlN2 + FgN2);
FtotN3 = (FlN3 + FgN3);
FtotN4 = (FlN4 + FgN4);
FtotN5 = (FlN5 + FgN5);

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),FtotN1)
% % plot(y(:,4),FtotN2)
% % plot(y(:,5),FtotN3)
% % plot(y(:,6),FtotN4)
% % plot(y(:,7),FtotN5)
% plot(y(:,3),-FtotN1,y(:,4),-FtotN2,y(:,5),-FtotN3,y(:,6),-FtotN4,y(:,7),-FtotN5)
% xlabel('Angle (degrees)')
% ylabel('Total Lateral Guide Forces (N)')
% title('Total Lateral Guide Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),-FtotN1)
% xlabel('Angle (degrees)')
% ylabel('Total Lateral Guide Force (N)')
% title('Total Lateral Guide Force (1 cylinder) - Angle')


% total tangential force
Ftott1 = Flt1 + Fgt1;
Ftott2 = Flt2 + Fgt2;
Ftott3 = Flt3 + Fgt3;
Ftott4 = Flt4 + Fgt4;
Ftott5 = Flt5 + Fgt5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Ftott1)
% % plot(y(:,4),Ftott2)
% % plot(y(:,5),Ftott3)
% % plot(y(:,6),Ftott4)
% % plot(y(:,7),Ftott5)
% plot(y(:,3),Ftott1,y(:,4),Ftott2,y(:,5),Ftott3,y(:,6),Ftott4,y(:,7),Ftott5)
% xlabel('Angle (degrees)')
% ylabel('Total Tangential Crankpin Forces (N)')
% title('Total Tangential Crankpin Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Ftott1)
% xlabel('Angle (degrees)')
% ylabel('Total Tangential Crankpin Force (N)')
% title('Total Tangential Crankpin Force (1 cylinder) - Angle')


% total radial force
FtotR1 = FlR1 + FgR1;
FtotR2 = FlR2 + FgR2;
FtotR3 = FlR3 + FgR3;
FtotR4 = FlR4 + FgR4;
FtotR5 = FlR5 + FgR5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),FtotR1)
% % plot(y(:,4),FtotR2)
% % plot(y(:,5),FtotR3)
% % plot(y(:,6),FtotR4)
% % plot(y(:,7),FtotR5)
% plot(y(:,3),FtotR1,y(:,4),FtotR2,y(:,5),FtotR3,y(:,6),FtotR4,y(:,7),FtotR5)
% xlabel('Angle (degrees)')
% ylabel('Total Radial Crankpin Forces (N)')
% title('Total Radial Crankpin Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),FtotR1)
% xlabel('Angle (degrees)')
% ylabel('Total Radial Force (N)')
% title('Total Radial Crankpin Force (1 cylinder) - Angle')


% total crankpin force
Ftots1 = Fls1 + Fgs1;
Ftots2 = Fls2 + Fgs2;
Ftots3 = Fls3 + Fgs3;
Ftots4 = Fls4 + Fgs4;
Ftots5 = Fls5 + Fgs5;

% figure 
% hold on
% subplot(1,2,1)
% % plot(y(:,3),Ftots1)
% % plot(y(:,4),Ftots2)
% % plot(y(:,5),Ftots3)
% % plot(y(:,6),Ftots4)
% % plot(y(:,7),Ftots5)
% plot(y(:,3),Ftots1,y(:,4),Ftots2,y(:,5),Ftots3,y(:,6),Ftots4,y(:,7),Ftots5)
% xlabel('Angle (degrees)')
% ylabel('Total Crankpin Forces (N)')
% title('Total Crankpin Forces - Angle')
% %figure
% subplot(1,2,2)
% plot(y(:,3),Ftots1)
% xlabel('Angle (degrees)')
% ylabel('Total Crankpin Force (N)')
% title('Total Crankpin Force (1 cylinder) - Angle')



% total Force
figure
hold on
plot(y(:,3),Fg1)
plot(y(:,3),Fl1)
plot(y(:,3),Ftot1)
legend({'Gas Force','Reciprocating Mass Inertia Force','Total Force'},'Location','northeast')
xlabel('Angle (degrees)')
ylabel('Piston Force (N)')
title('Piston Force - Angle')


% total lateral guide Force
figure
hold on
plot(y(:,3),-FgN1)
plot(y(:,3),-FlN1)
plot(y(:,3),-FtotN1)
legend({'Gas Lateral Guide Force','Reciprocating Mass Inertia Lateral Guide Force','Total Lateral Guide Force'},'Location','northeast')
xlabel('Angle (degrees)')
ylabel('Lateral Guide Force')
title('Lateral Guide Force - Angle')


% total crankpin  Force
figure
hold on
plot(y(:,3),Fgs1)
plot(y(:,3),Fls1)
plot(y(:,3),Ftots1)
legend({'Gas Crankpin Force','Reciprocating Mass Inertia Crankpin Force','Total Crankpin Force'},'Location','northeast')
xlabel('Angle (degrees)')
ylabel('Crankpin Force (N)')
title('Crankpin Force - Angle')


% total tangential  Force
figure
hold on
plot(y(:,3),Fgt1)
plot(y(:,3),Flt1)
plot(y(:,3),Ftott1)
legend({'Gas Tangential Crankpin Force','Reciprocating Mass Inertia Tangential Crankpin Force','Total Tangential Crankpin Force'},'Location','northeast')
xlabel('Angle (degrees)')
ylabel('Tangential Crankpin Force (N)')
title('Tangential Crankpin Force - Angle')


% total radial  Force
figure
hold on
plot(y(:,3),FgR1)
plot(y(:,3),FlR1)
plot(y(:,3),FtotR1)
legend({'Gas Radial Crankpin Force','Reciprocating Mass Inertia Radial Crankpin Force','Total Radial Crankpin Force'},'Location','northeast')
xlabel('Angle (degrees)')
ylabel('Radial Crankpin Force (N)')
title('Radial Crankpin Force - Angle')


























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Barred Speed Range
% lamda pou temnoun tis idiosixnotites
% idiosixnotika 1 = 330.08, lamda: 4 5 6 7
% idiosixnotika 2 = 502.10, lamda: 6 7 8 9 10 11 12 13 14 15 16 17 18
% idiosixnotika 3 = 1872.20, lamda: 20 21 22 23 24 25
% pano apo tin idiosixnotita 3 den temnetai me kapoio lamda


% figure
% hold on
% 
% n = [0:110];
% nf1 = idiosixnotites_rpm2(2)*ones(size(n));
% nf2 = idiosixnotites_rpm2(3)*ones(size(n));
% nf3 = idiosixnotites_rpm2(4)*ones(size(n));
% 
% plot(n,nf1)
% plot(n,nf2)
% plot(n,nf3)
% 
% lamda1 = 1;
% lamda2 = 2;
% lamda3 = 3;
% lamda4 = 4;
% lamda5 = 5;
% lamda6 = 6;
% lamda7 = 7;
% lamda8 = 8;
% lamda9 = 9;
% lamda10 = 10;
% lamda11 = 11;
% lamda12 = 12;
% lamda13 = 13;
% lamda14 = 14;
% lamda15 = 15;
% lamda16 = 16;
% lamda17 = 17;
% lamda18 = 18;
% lamda19 = 19;
% lamda20 = 20;
% lamda21 = 21;
% lamda22 = 22;
% lamda23 = 23;
% lamda24 = 24;
% 
% l1 = lamda1*n;
% l2 = lamda2*n;
% l3 = lamda3*n;
% l4 = lamda4*n;
% l5 = lamda5*n;
% l6 = lamda6*n;
% l7 = lamda7*n;
% l8 = lamda8*n;
% l9 = lamda9*n;
% l10 = lamda10*n;
% l11 = lamda11*n;
% l12 = lamda12*n;
% l13 = lamda13*n;
% l14 = lamda14*n;
% l15 = lamda15*n;
% l16 = lamda16*n;
% l17 = lamda17*n;
% l18 = lamda18*n;
% l19 = lamda19*n;
% l20 = lamda20*n;
% l21 = lamda21*n;
% l22 = lamda22*n;
% l23 = lamda23*n;
% l24 = lamda24*n;
% 
% plot(n,l1)
% plot(n,l2)
% plot(n,l3)
% plot(n,l4)
% plot(n,l5)
% plot(n,l6)
% plot(n,l7)
% plot(n,l8)
% plot(n,l9)
% plot(n,l10)
% plot(n,l11)
% plot(n,l12)
% plot(n,l13)
% plot(n,l14)
% plot(n,l15)
% plot(n,l16)
% plot(n,l17)
% plot(n,l18)
% plot(n,l19)
% plot(n,l20)
% plot(n,l21)
% plot(n,l22)
% plot(n,l23)
% plot(n,l24)
% 
% xlabel('Angular Velocity (rpm)')
% ylabel('Natural Frequencies (rpm)')
% title('Campbell Diagram')
% 
% idio_nkr = idiosixnotites_rpm2(2)/lamda5; % krataw to lamda=5, z(cylinder number)=5, kratao tin proti idiosixnotita
% idio_rkr = idio_nkr/NMCR;
% 
% % apagoreumeni_perioxi_strofon (a.p.s.)
% aps_min = 16*idio_nkr/(18-idio_rkr);
% aps_max = (18-idio_rkr)*idio_nkr/16;
% 
% 
% 
% intemidiate_shaft_torque = K(10,11)*(y(:,11)-y(:,10));
% tail_shaft_torque =  K(11,12)*(y(:,12)-y(:,11));
% 
% smooth_intemidiate_shaft_torque = smoothdata(smoothdata(smoothdata(smoothdata(intemidiate_shaft_torque))));
% smooth_tail_shaft_torque = smoothdata(smoothdata(smoothdata(smoothdata(tail_shaft_torque))));
% 
% wkr = 5.655;
% nkr = idiosixnotites_rpm2(2)/5;
% %wkr = nkr*2*pi/60;
% 
% wintsh = y22_smooth;
% wtailsh = y23_smooth;
% 
% QR = 7/100;
% 
% Qintsh = 1./sqrt((1-(wintsh/wkr).^2).^2 + (2*QR)*(wintsh/wkr).^2);
% Qtailsh = 1./sqrt((1-(wtailsh/wkr).^2).^2 + (2*QR)*(wtailsh/wkr).^2);
% 
% uts_intsh = 800; % N/mm^2
% uts_tailsh = 600; % N/mm^2
% 
% do_intsh = 340; % mm
% do_tailsh = 450; % mm
% 
% ck_intsh = 1;
% ck_tailsh = 0.55;
% 
% cD_intsh = 0.35+0.93*do_intsh^(-0.2);
% cD_tailsh = 0.35+0.93*do_tailsh^(-0.2);
% 
% nmcr=99;
% wmcr = 2*pi*nmcr/60;
% 
% w_start = start_rpm*2*pi/60;
% 
% wintsh1 = w_start:0.1:0.9*wmcr;
% wintsh2 = 0.9*wmcr:0.1:1.05*wmcr;
% 
% tc_intsh_1 = (uts_intsh+160)/18*ck_intsh*cD_intsh*(3-2*(wintsh1/wmcr).^2);
% tc_intsh_2 = ones(size(wintsh2))*(uts_intsh+160)/18*ck_intsh*cD_intsh*1.38;
% tt_intsh_1 = 1.7*tc_intsh_1/sqrt(ck_intsh);
% tt_intsh_2 = 1.7*tc_intsh_2/sqrt(ck_intsh);
% 
% wtailsh1 = w_start:0.1:0.9*wmcr;
% wtailsh2 = 0.9*wmcr:0.1:1.05*wmcr;
% 
% tc_tailsh_1 = (uts_tailsh+160)/18*ck_tailsh*cD_tailsh*(3-2*(wtailsh1/wmcr).^2);
% tc_tailsh_2 = ones(size(wtailsh2))*(uts_tailsh+160)/18*ck_tailsh*cD_tailsh*1.38;
% tt_tailsh_1 = 1.7*tc_tailsh_1/sqrt(ck_tailsh);
% tt_tailsh_2 = 1.7*tc_tailsh_2/sqrt(ck_tailsh);
% 
% static_stress_intsh = (smooth_intemidiate_shaft_torque*((do_intsh/1000)/2)*32/(pi*(do_intsh/1000)^4))/10^6;
% static_stress_tailsh = (smooth_tail_shaft_torque*((do_tailsh/1000)/2)*32/(pi*(do_tailsh/1000)^4))/10^6;
% 
% 
% dynamic_stress_intsh = static_stress_intsh.*Qintsh;
% dynamic_stress_tailsh = static_stress_tailsh.*Qtailsh;
% 
% 
% figure
% hold on
% h(1) = plot(60*wintsh1/(2*pi),tc_intsh_1,'k')
% h(2) = plot(60*wintsh2/(2*pi),tc_intsh_2,'k')
% h(3) = plot(60*wintsh1/(2*pi),tt_intsh_1,'r')
% h(4) = plot(60*wintsh2/(2*pi),tt_intsh_2,'r')
% h(5) = plot(60*wintsh/(2*pi),dynamic_stress_intsh,'b')
% %plot(60*wintsh/(2*pi),static_stress_intsh,'y')
% legend(h([1 3 5]),'Continuous Operation Limit','Transient Operation Limit','Actual Stress','Location','northeast')
% 
% xlabel('Angle Velocity (rpm)')
% ylabel('Torsional Stress (N/mm^2)')
% title('Intermidiate Shaft Torsional Stress')
% 
% 
% figure
% hold on
% hh(1) = plot(60*wtailsh1/(2*pi),tc_tailsh_1,'k')
% hh(2) = plot(60*wtailsh2/(2*pi),tc_tailsh_2,'k')
% hh(3) = plot(60*wtailsh1/(2*pi),tt_tailsh_1,'r')
% hh(4) = plot(60*wtailsh2/(2*pi),tt_tailsh_2,'r')
% hh(5) = plot(60*wtailsh/(2*pi),dynamic_stress_tailsh,'b')
% %plot(60*wtailsh/(2*pi),static_stress_tailsh,'y')
% legend(hh([1 3 5]),'Continuous Operation Limit','Transient Operation Limit','Actual Stress','Location','northeast')
% xlabel('Angle Velocity (rpm)')
% ylabel('Torsional Stress (N/mm^2)')
% title('Tail Shaft Torsional Stress')







