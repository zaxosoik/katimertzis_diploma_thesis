function p_interpol = p_interpol(y,load)
cyl = [3 6 5 4 7];  % seira kaysis kilindron DoF

P = zeros(12,1);  % arxikopoiisi

% diafora fasis kilindron
n = 5;
aa = 360/n;
a = aa*pi/180;
a = [0 1*a 2*a 3*a 4*a];
if load==25
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
elseif load ==50
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
elseif load ==75
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

else
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

end
p_interpol=P;
end