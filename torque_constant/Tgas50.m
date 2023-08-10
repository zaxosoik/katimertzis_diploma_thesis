function Tgas50 = Tgas50(y,Ds)  % ropi logo kaysis

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
            P1 = 154888.88*thita^3 + 1101111.97*thita^2 + 2598942.5*thita + 2281157.77;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 3007037.78*thita^3 + 11515997.52*thita^2 + 15449184.93*thita + 7657601.24;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -149755092.04*thita^5 - 373681467.88*thita^4 - 347745253.65*thita^3 - 131831332.17*thita^2 - 2572554.98*thita + 9534248.33;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = -521291354.75*thita^3 + 69916940.96*thita^2 + 19561879.7*thita + 10069151.9;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = -151819660.78*thita^6 + 701539476.81*thita^5 - 1319093514.04*thita^4 + 1272750060.74*thita^3 - 637308115.41*thita^2 + 132366862.88*thita + 3501064.82;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = -1245949.93*thita^3 + 6900191.48*thita^2 - 13268462.86*thita + 9524157.28;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = -1428125.23*thita^4 + 14408298.52*thita^3 - 53496237.02*thita^2 + 86117734.26*thita - 50028400.3; 
        end
        
    else
        thita = WW(ii)-pi;
        if  (-180*2*pi/360 <= thita) && (thita <= -75*2*pi/360)
            P1 = 154888.88*thita^3 + 1101111.97*thita^2 + 2598942.5*thita + 2281157.77;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 3007037.78*thita^3 + 11515997.52*thita^2 + 15449184.93*thita + 7657601.24;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -149755092.04*thita^5 - 373681467.88*thita^4 - 347745253.65*thita^3 - 131831332.17*thita^2 - 2572554.98*thita + 9534248.33;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = -521291354.75*thita^3 + 69916940.96*thita^2 + 19561879.7*thita + 10069151.9;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = -151819660.78*thita^6 + 701539476.81*thita^5 - 1319093514.04*thita^4 + 1272750060.74*thita^3 - 637308115.41*thita^2 + 132366862.88*thita + 3501064.82;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = -1245949.93*thita^3 + 6900191.48*thita^2 - 13268462.86*thita + 9524157.28;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = -1428125.23*thita^4 + 14408298.52*thita^3 - 53496237.02*thita^2 + 86117734.26*thita - 50028400.3; 
        end
    end
    
       P(cyl(ii)) = P1;

end


B = 0.5;     % bore (m)
l = 2.05;   % piston rod
r = 2.05/2;  % crankshaft, stroke/2

%Tgas50 = - 0.5*(78.6/62.4)*(3175/2550)*Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
%Tgas50=P;  (mono gia plot kateutheian tin piesi)
Tgas50 = - Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")

end




% function Tgas50 = Tgas50(y,Ds)  % ropi logo kaysis
% 
% % sintelestes polionimou eksisosis piesis me gonia, ypologismeno apo excel
% % proseggisi kampilis piesis me 2 polionima
% A1 = -387111.271;
% A2 = -3516443.323;
% A3 = -11242532.833;
% A4 = -12853241.478;
% A5 = 4106532.012;
% A6 = 18508384.141;
% A7 = 10184448.199;
% 
% B1 = -396078.051;
% B2 = 3733101.994;
% B3 = -12245582.647;
% B4 = 13254354.028;
% B5 = 10241683.203;
% B6 = -31769735.321;
% B7 = 19007153.862;
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
% Tgas50 = - Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
% %Tgas50=P;  (mono gia plot kateutheian tin piesi)
% end
% 
% % to "-" sto Ds einai logo tis gonias apo 0 ews 360
% % gia gonia apo -180 ews 180 den xreiazetai "-"