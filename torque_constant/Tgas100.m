function Tgas100 = Tgas100(y,Ds)  % ropi logo kaysis

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
        if  (-180*2*pi/360 <= thita) && (thita<= -75*2*pi/360)
            P1 = 284505.848*thita.^4 + 2696825.148*thita.^3 + 9390615.502*thita.^2 + 14312416.462*thita + 8546183.518;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 4855647.451*thita.^3 + 18575778.63*thita.^2 + 24899462.397*thita + 12336199.891;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -264604096.963*thita.^5 - 648141723.985*thita.^4 - 594338908.642*thita.^3 - 224278706.884*thita.^2 - 6281780.506*thita + 15149555.011;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = 21567932251.5*thita.^5 - 9533496708.344*thita.^4 + 1171491680.707*thita.^3 - 13928528.891*thita.^2 + 1939307.466*thita + 15284529.245;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = 61210688.219*thita.^6 - 156028953.611*thita.^5 + 42977050.473*thita.^4 + 209095009.908*thita.^3 - 227428427.165*thita.^2 + 62500954.834*thita + 10934574.109;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = -2094460.509*thita.^3 + 11604432.739*thita.^2 - 22328360.235*thita + 16047867.726;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = -1180785.609*thita.^4 + 13097773.143*thita.^3 - 53291633.829*thita.^2 + 93617836.383*thita - 58852082.548; 
        end
        
    else
        thita = WW(ii)-pi;
        if  (-180*2*pi/360 <= thita) && (thita<= -75*2*pi/360)
            P1 = 284505.848*thita.^4 + 2696825.148*thita.^3 + 9390615.502*thita.^2 + 14312416.462*thita + 8546183.518;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 4855647.451*thita.^3 + 18575778.63*thita.^2 + 24899462.397*thita + 12336199.891;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -264604096.963*thita.^5 - 648141723.985*thita.^4 - 594338908.642*thita.^3 - 224278706.884*thita.^2 - 6281780.506*thita + 15149555.011;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = 21567932251.5*thita.^5 - 9533496708.344*thita.^4 + 1171491680.707*thita.^3 - 13928528.891*thita.^2 + 1939307.466*thita + 15284529.245;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = 61210688.219*thita.^6 - 156028953.611*thita.^5 + 42977050.473*thita.^4 + 209095009.908*thita.^3 - 227428427.165*thita.^2 + 62500954.834*thita + 10934574.109;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = -2094460.509*thita.^3 + 11604432.739*thita.^2 - 22328360.235*thita + 16047867.726;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = -1180785.609*thita.^4 + 13097773.143*thita.^3 - 53291633.829*thita.^2 + 93617836.383*thita - 58852082.548; 
        end
    end
    
       P(cyl(ii)) = P1;

end


B = 0.5;     % bore (m)
l = 2.05;   % piston rod
r = 2.05/2;  % crankshaft, stroke/2

%Tgas=P;  (mono gia plot kateutheian tin piesi)
Tgas100 = - Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")

end




% function Tgas = Tgas(y,Ds)  % ropi logo kaysis
% 
% % sintelestes polionimou eksisosis piesis me gonia, ypologismeno apo excel
% % proseggisi kampilis piesis me 2 polionima
% A1 = -1097738.273;
% A2 = -10471290.719;
% A3 = -36895735.523;
% A4 = -56066717.468;
% A5 = -25933574.076;
% A6 = 16795306.331;
% A7 = 14855855.738;
% 
% B1 = -903155.971;
% B2 = 9458772.977;
% B3 = -37598904.652;
% B4 = 68248722.917;
% B5 = -48341409.942;
% B6 = -7029422.81;
% B7 = 19350434.208;
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
% Tgas = - Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
% %Tgas=P;  (mono gia plot kateutheian tin piesi)
% end
% 
% % to "-" sto Ds einai logo tis gonias apo 0 ews 360
% % gia gonia apo -180 ews 180 den xreiazetai "-"
% 
% 
