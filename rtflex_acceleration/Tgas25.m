function Tgas25 = Tgas25(y,Ds)  % ropi logo kaysis

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
            P1 = 80313.61*thita^3 + 625010.23*thita^2 + 1627039.23*thita + 1572902.76;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 2383875.73*thita^3 + 9110983.31*thita^2 + 12202094.95*thita + 6042057.33;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -72030806.37*thita^4 - 120615576.24*thita^3 - 57033032.47*thita^2 + 3872426*thita + 7689975.82;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = -309878347.71*thita^3 + 17772401.98*thita^2 + 15566699.56*thita + 7996142.19;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = 99020332.96*thita^5 - 344074367.19*thita^4 + 451073482.96*thita^3 - 262860468.86*thita^2 + 51616469.01*thita + 6559543.19;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = 2345992.28*thita^3 - 9091048.2*thita^2 + 9832747.41*thita - 1794887.14;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = -169816.37*thita^3 + 1329584.4*thita^2 - 3407610.6*thita + 2989843.33; 
        end
        
    else
        thita = WW(ii)-pi;
        if  (-180*2*pi/360 <= thita) && (thita <= -75*2*pi/360)
            P1 = 80313.61*thita^3 + 625010.23*thita^2 + 1627039.23*thita + 1572902.76;
        elseif (-75*2*pi/360 < thita) && (thita <= -40*2*pi/360)
            P1 = 2383875.73*thita^3 + 9110983.31*thita^2 + 12202094.95*thita + 6042057.33;
        elseif (-40*2*pi/360 < thita) && (thita <= -2*2*pi/360)
            P1 = -72030806.37*thita^4 - 120615576.24*thita^3 - 57033032.47*thita^2 + 3872426*thita + 7689975.82;
        elseif (-2*2*pi/360 < thita) && (thita <= 10*2*pi/360)
            P1 = -309878347.71*thita^3 + 17772401.98*thita^2 + 15566699.56*thita + 7996142.19;
        elseif (10*2*pi/360 < thita) && (thita <= 60*2*pi/360)
            P1 = 99020332.96*thita^5 - 344074367.19*thita^4 + 451073482.96*thita^3 - 262860468.86*thita^2 + 51616469.01*thita + 6559543.19;
        elseif (60*2*pi/360 < thita) && (thita <= 110*2*pi/360)
            P1 = 2345992.28*thita^3 - 9091048.2*thita^2 + 9832747.41*thita - 1794887.14;
        elseif (110*2*pi/360 < thita) && (thita <= 180*2*pi/360)
            P1 = -169816.37*thita^3 + 1329584.4*thita^2 - 3407610.6*thita + 2989843.33; 
        end
    end
    
       P(cyl(ii)) = P1;

end


B = 0.5;     % bore (m)
l = 2.05;   % piston rod
r = 2.05/2;  % crankshaft, stroke/2

Tgas25 = -Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
%Tgas25=P;  (mono gia plot kateutheian tin piesi)
end




% function Tgas25 = Tgas25(y,Ds)  % ropi logo kaysis
% 
% % sintelestes polionimou eksisosis piesis me gonia, ypologismeno apo excel
% % proseggisi kampilis piesis me 2 polionima
% A1 = -318938.831;
% A2 = -2908095.157;
% A3 = -9369554.519;
% A4 = -11030778.528;
% A5 = 2433160.508;
% A6 = 14231568.299;
% A7 = 7959816.042;
% 
% B1 = -129116.270;
% B2 = 974436.809;
% B3 = -1643547.869;
% B4 = -4697241.121;
% B5 = 21025998.062;
% B6 = -28755495.014;
% B7 = 14512126.484;
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
% Tgas25 = - Ds*(pi/4)*(B^2)*P;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")
% %Tgas25=P;  (mono gia plot kateutheian tin piesi)
% end
% 
% % to "-" sto Ds einai logo tis gonias apo 0 ews 360
% % gia gonia apo -180 ews 180 den xreiazetai "-"


