% 1os tropos, excel friction-damping, ypoektimisi 11%
% function Tfriction = Tfriction(wmega)  % ropi logo trivwn
% 
% Stroke = 2.05;
% Z = 5; % cylinder number
% D = 0.5;
% VH = Z*Stroke*pi*D^2/4;
% C1 = 5258.3;
% F1 = - VH * C1 / (2*pi*Z);
% Tfr1 = [0 0 F1 F1 F1 F1 F1 0 0 0 0 0]';
% 
% C2 = 7779.4;
% F2 = VH*C2 / (2*pi*Z);
% Tfr222 = [0 0 F2 F2 F2 F2 F2 0 0 0 0 0];
% Tfr22 = diag(Tfr222);
% Tfr2 = Tfr22*wmega;
% 
% Tfriction = Tfr1 + Tfr2;
% 
% end



% % 2os tropos, pfr = 0.137 + 0.005*pimax + 0.162*c, yperektimisi 36%
% function Tfriction = Tfriction(wmega)  % ropi logo trivwn
% 
% Stroke = 2.05;
% Z = 5; % cylinder number
% D = 0.5;
% VH = Z*Stroke*pi*D^2/4;
% pimax = 19.97; % maximum mean indicated pressure (bar)
% F1 = ((VH/(2*pi))*(0.137+0.005*pimax)*(10^5))/Z;
% Tfr1 = [0 0 F1 F1 F1 F1 F1 0 0 0 0 0]';
% 
% 
% F2 = (VH/(2*pi))*0.162*Stroke*(10^5)/(Z*pi);
% Tfr222 = [0 0 F2 F2 F2 F2 F2 0 0 0 0 0];
% Tfr22 = diag(Tfr222);
% Tfr2 = Tfr22*wmega;
% 
% Tfriction = Tfr1 + Tfr2;
% 
% end



% 3os tropos, absolute damping coefficient
function Tfriction = Tfriction(wmega)  % ropi logo trivwn

% sintelestes gia Tfriction - Absolute Damping
C1 = 540;
C2 = 0;
C = [0 0 C1 C1 C1 C1 C1 0 0 0 0 C2];

Cf = diag(C);  % Cf diagonios pinakas 16x16

Tfriction = Cf*wmega;

end



% 4os tropos, Pfr = A*n^b
% function Tfriction = Tfriction(wmega)  % ropi logo trivwn
% 
% % sintelestes gia Pfr = A*n^b
% b = 1; % b = 1-1.2 (megales mixanes:1, mikres mixanes:1.2)
% 
% pfr = 0.85*10^5; %pfriction MCR : 0.85 bar
% Stroke = 2.05;
% Z = 5; % cylinder number
% D = 0.5;
% VH = Z*Stroke*pi*D^2/4;
% NMCR = 99;
% 
% A = (pfr*VH*NMCR/60)/(NMCR^b);
% 
% AA = A*(60/(2*pi))^b;
% AAA = AA/Z;
% 
% C = [0 0 AAA AAA AAA AAA AAA 0 0 0 0 0];
% 
% Cf = diag(C);  % Cf diagonios pinakas 16x16
% 
% Tfriction = Cf*(60*wmega/2*pi).^(b-1);
% end



% % 5os tropos, fmep, yperektimisi 26%
% function Tfriction = Tfriction(wmega)  % ropi logo trivwn
% 
% Stroke = 2.05;
% Z = 5; % cylinder number
% D = 0.5;
% VH = Z*Stroke*pi*D^2/4;
% 
% F1 = ((VH/(2*pi))*(0.0384*(1+1/Z)))*10^6/Z;
% Tfr1 = [0 0 F1 F1 F1 F1 F1 0 0 0 0 0]';
% 
% F2 = ((VH/(2*pi))*(3/(D*1000)))*10^6/Z;
% Tfr2 = [0 0 F2 F2 F2 F2 F2 0 0 0 0 0]';
% 
% F3A = (VH/(2*pi))*0.018*0.0179*10^6/Z;
% Tfr3A3A3A = [0 0 F3A F3A F3A F3A F3A 0 0 0 0 0];
% Tfr3A3A = diag(Tfr3A3A3A);
% Tfr3A = Tfr3A3A*wmega.^2;
% 
% F3B = -(VH/(2*pi))*0.018*0.0005*10^6/Z;
% Tfr3B3B3B = [0 0 F3B F3B F3B F3B F3B 0 0 0 0 0];
% Tfr3B3B = diag(Tfr3B3B3B);
% Tfr3B = Tfr3B3B*wmega;
% 
% F3C = -((VH/(2*pi))*0.018*4*10^(-5))*10^6/Z;
% Tfr3C = [0 0 F3C F3C F3C F3C F3C 0 0 0 0 0]';
% 
% Tfr3 = Tfr3A + Tfr3B + Tfr3C;
% 
% F4 = (VH/(2*pi))*0.004*(Stroke/pi)*10^6/Z;
% Tfr444 = [0 0 F4 F4 F4 F4 F4 0 0 0 0 0];
% Tfr44 = diag(Tfr444);
% Tfr4 = Tfr44*wmega;
% 
% Tfriction = Tfr1 + Tfr2 + Tfr3 + Tfr4;
% 
% end