function Tgas = Tgas(y,Ds,time,t_metabolis,load1,load2)  % ropi logo kaysis
if time<t_metabolis
    if load1==25
        Tgas=Tgas25(y,Ds);
    elseif load1 ==50
        Tgas=Tgas50(y,Ds);
    elseif load1 ==75
        Tgas=Tgas75(y,Ds);
    else
        Tgas=Tgas100(y,Ds);
    end
else
    if load2==25
        Tgas=Tgas25(y,Ds);
    elseif load2 ==50
        Tgas=Tgas50(y,Ds);
    elseif load2 ==75
        Tgas=Tgas75(y,Ds);
    else
        Tgas=Tgas100(y,Ds);
    end
end

end

