function Tgas = Tgas(y,Ds,rpm)  % ropi logo kaysis

B = 0.5;     % bore (m)

%load = [25,50,75,100];
load = [25,50,75];


%w = 2*pi*[59.3,79.4,89.7,98.9]/60; % RPM gia 25,50,75,100% load
w = 2*pi*[59.3,79.4,89.7]/60; % RPM gia 25,50,75,100% load

% orizoume to p_load_final, dhladh tin teliki eksisosi piesis 
for i = 1:length(w)
    p_load(:,i) = p_interpol(y,load(i));
end

p_load_final = p_load(:,1); % orizoume genika tin piesi isi me 25% fortio
                            % sto apo kato forloop pairnoume periptoseis
                            % den meletame dhladh katholou me grammiki paremvoli gia strofes mikroteres tou 25%
                            % gia mikroteres strofes theoroume fortio 25%


for i = 1:length(w)-1 %w(1), w(2), w(3), 
    
    if rpm >= w(i) && rpm < w(i+1) % praktika an einai anamesa sto w(1)=62.4 kai sto w(3+1=4)=99 kanei grammiki paramvoli
        
        p_load_final = perithorio(rpm,2*pi*98.9/60)*((p_load(:,i)+(-w(i)+rpm)/(w(i+1)-w(i))*(p_load(:,i+1)-p_load(:,i))));
        
        break
    
    elseif i == length(w)-1 && rpm >= w(i+1) % an einai megalitero apo  w(3+1=4)=99 krataw tin timi gia tin piesi gia 99
                                             % edo mporo na epilekso se poio fortio na stamatisei
                                             % allazontas sto rpm >= w(i+1) vazontaw poio w thelo (p.x. w(i), w(i-1) ...)
                                             % allazo episis to teliko fortio p_load(:,4) apo kato, (allazo to ''4'')

        
       % p_load_final = p_load(:,3);
              p_load_final = p_load(:,length(w));

        
        break
    
    end
end

Tgas = - Ds*(pi/4)*(B^2)*p_load_final;  % eksisosi ropis logo kaysis (epeksigisi sto main gia to "-")

end

% mporoume na ksekinisoume to sistima apo opoio rpm theloume megalitero tou
% (rpm=25% load) (opoiadipote tyxaia timi)
% stamatame to programma se sigkekrimeno fortio (25%,50%,75%,100%) analoga
% ti tou orisoume sto elseif sto deutero forloop
