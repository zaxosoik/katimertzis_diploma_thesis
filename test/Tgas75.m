function Tgas75 = Tgas75(y,Ds)  % ropi logo kaysis

cyl = [3 6 5 4 7];  % seira kaysis kilindron DoF

P = zeros(12,1);  % arxikopoiisi

% diafora fasis kilindron
n = 5;
aa = 360/n;
a = aa*pi/180;
a = [0 1*a 2*a 3*a 4*a];


for ii = 1:length(cyl)
    
    WW(ii) = y(cyl(ii)) + a(ii);
    
    if WW(ii) > 2*pi
        thita = mod(WW(ii),2*pi)-pi;   % metaferei tis gonies me vasi to diagramma Cylinder-Pressure se -180 eos 180 giati etsi kataskeuastike to diagramma Cylinder Pressure (deinamodeiktiko) sto excel
        if  (-180*2*pi/360 <= thita) && (thita <= -75*2*pi/360)
            P1 = 207349.4*thita^3 + 1465230.91*thita^2 + 3437339.32*thita + 3022724.75;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 4114632.36*thita^3 + 15752557.36*thita^2 + 21127840.98*thita + 10472098.15;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -227261453.88*thita^5 - 555865560.96*thita^4 - 508895631.33*thita^3 - 191656854.83*thita^2 - 5295599.06*thita + 12927359.97;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = -499210569.83*thita^3 + 96577737.61*thita^2 + 11881604.15*thita + 13191446.26;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = -76879972.11*thita^6 + 430748250.34*thita^5 - 952205781.03*thita^4 + 1053923180.69*thita^3 - 591628720.85*thita^2 + 133092224.47*thita + 5530832.93;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = -1664464.52*thita^3 + 9219481.67*thita^2 - 17730997.22*thita + 12729039.79;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = 553098.37*thita^3 - 3657787.96*thita^2 + 7178784.76*thita - 3279913.79; 
        end
        
    else
        thita = WW(ii)-pi;
        if  (-180*2*pi/360 <= thita) && (thita <= -75*2*pi/360)
            P1 = 207349.4*thita^3 + 1465230.91*thita^2 + 3437339.32*thita + 3022724.75;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 4114632.36*thita^3 + 15752557.36*thita^2 + 21127840.98*thita + 10472098.15;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -227261453.88*thita^5 - 555865560.96*thita^4 - 508895631.33*thita^3 - 191656854.83*thita^2 - 5295599.06*thita + 12927359.97;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = -499210569.83*thita^3 + 96577737.61*thita^2 + 11881604.15*thita + 13191446.26;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = -76879972.11*thita^6 + 430748250.34*thita^5 - 952205781.03*thita^4 + 1053923180.69*thita^3 - 591628720.85*thita^2 + 133092224.47*thita + 5530832.93;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = -1664464.52*thita^3 + 9219481.67*thita^2 - 17730997.22*thita + 12729039.79;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = 553098.37*thita^3 - 3657787.96*thita^2 + 7178784.76*thita - 3279913.79; 
        end
    end
    
       P(cyl(ii)) = P1;

end


B = 0.5;     % bore (m)
l = 2.05;   % piston rod
r = 2.05/2;  % crankshaft, stroke/2

Tgas75 = -Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
%Tgas75=P;  (mono gia plot kateutheian tin piesi)
end




% function Tgas75 = Tgas75(y,Ds)  % ropi logo kaysis
% 
% % sintelestes polionimou eksisosis piesis me gonia, ypologismeno apo excel
% % proseggisi kampilis piesis me 2 polionima
% A1 = -742273.181;
% A2 = -6965427.022;
% A3 = -23808384.909;
% A4 = -33506006.09;
% A5 = -9164963.05;
% A6 = 19304175.851;
% A7 = 13191700.648;
% 
% B1 = -663101.499;
% B2 = 6613317.576;
% B3 = -24243012.002;
% B4 = 37038518.888;
% B5 = -11277146.419;
% B6 = -26558094.212;
% B7 = 21544361.95;
% 
% 
% cyl = [3 6 5 4 7];  % seira kaysis kilindron DoF
% 
% max_pressure_angle_degrees = 10;
% max_pressure_angle = max_pressure_angle_degrees*2*pi/360;  % rad
% 
% % ypologismos eksisosis piesis
% P = zeros(12,1);  % arxikopoiisi
% 
% % diafora fasis kilindron
% n = 5;
% aa = 360/n;
% a = aa*pi/180;
% 
% a = [0 1*a 2*a 3*a 4*a];
% 
% for ii = 1:length(cyl)
%     
%     WW(ii) = y(cyl(ii)) + a(ii);
% 
%     if WW(ii) > 2*pi
%         thita = mod(WW(ii),2*pi)-pi;   % metaferei tis gonies me vasi to diagramma Cylinder-Pressure se -180 eos 180 giati etsi kataskeuastike to diagramma Cylinder Pressure (deinamodeiktiko) sto excel
%         if thita > max_pressure_angle  %gonia megistis piesis sto diagramma, allagi eksisosis sto excel
%             P1 = B1*thita^6+B2*thita^5+B3*thita^4+B4*thita^3+B5*thita^2+B6*thita+B7;
%         else
%             P1 = A1*thita^6+A2*thita^5+A3*thita^4+A4*thita^3+A5*thita^2+A6*thita+A7;
%         end
%     else
%         thita = WW(ii)-pi;
%         if thita > max_pressure_angle
%             P1 = B1*thita^6+B2*thita^5+B3*thita^4+B4*thita^3+B5*thita^2+B6*thita+B7;
%         else
%             P1 = A1*thita^6+A2*thita^5+A3*thita^4+A4*thita^3+A5*thita^2+A6*thita+A7;
%         end
%     end
%    P(cyl(ii)) = P1;
% 
% end
% 
% 
% B = 0.5;     % bore (m)
% l = 2.05;   % piston rod
% r = 2.05/2;  % crankshaft, stroke/2
% 
% Tgas75 = - Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
% %Tgas75=P;  (mono gia plot kateutheian tin piesi)
% end
% 
% % to "-" sto Ds einai logo tis gonias apo 0 ews 360
% % gia gonia apo -180 ews 180 den xreiazetai "-"


