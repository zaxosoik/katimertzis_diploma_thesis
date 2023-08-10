function Tgas = Tgas(y,Ds,rpm)  % ropi logo kaysis
B = 0.5;     % bore (m)

load =[25,50,75,100];
w =2*pi*[59.37,79.7,89.96,99.08]/60; % RPM gia 25,50,75,100% load
%w = 2*pi*w/60;
for i=1:length(w)
    p_load(:,i) = p_interpol(y,load(i));
end
 p_load_final = p_load(:,1);
for i=1:length(w)-1
    if rpm>=w(i) && rpm<w(i+1) 
        %if rpm>=96*2*pi/60
        %    p_load_final = p_load(:,4);
        %    break
        %else
            %p_load_final = perithorio(rpm,2*pi*99.08/60)*(((w(i+1)-rpm)/(w(i+1)-w(i))*p_load(:,i)+(-w(i)+rpm)/(w(i+1)-w(i))*p_load(:,i+1)));
            p_load_final = perithorio(rpm,2*pi*99.08/60)*((p_load(:,i)+(-w(i)+rpm)/(w(i+1)-w(i))*(p_load(:,i+1)-p_load(:,i))));

            %p_load_final = 1*(((w(i+1)-rpm)/(w(i+1)-w(i))*p_load(:,i)+(-w(i)+rpm)/(w(i+1)-w(i))*p_load(:,i+1)));
            break
        %end
    %elseif i==1 && rpm<w(1)
    %    p_load_final =p_load(:,1);
    %    break
    
    elseif i==length(w)-1 && rpm>=w(i+1)
        p_load_final = p_load(:,4);
        break
    end
end
%p_load_final = perithorio(rpm,6350000,99*2*pi/60)*p_load_final;
%p_load_final = 2*p_load_final;
%Tgas=P;  (mono gia plot kateutheian tin piesi)
Tgas = - Ds*(pi/4)*(B^2)*p_load_final;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")

%if time<t_metabolis
%    if load1==25
%        Tgas=Tgas25(y,Ds);
%    elseif load1 ==50
%        Tgas=Tgas50(y,Ds);
%    elseif load1 ==75
%        Tgas=Tgas75(y,Ds);
%    else
%        Tgas=Tgas100(y,Ds);
%    end
%else
%    if load2==25
%        Tgas=Tgas25(y,Ds);
%    elseif load2 ==50
%        Tgas=Tgas50(y,Ds);
%    elseif load2 ==75
%        Tgas=Tgas75(y,Ds);
%    else
%        Tgas=Tgas100(y,Ds);
%    end
%end

end

