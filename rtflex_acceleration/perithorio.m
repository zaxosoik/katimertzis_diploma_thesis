% se auto to .m apla ftiaxnoume to perithorio isxuos

function perithorio = perithorio(rpm,n_MCR)

if rpm > 59*2*pi/60 && rpm < 95*2*pi/60 % oria rpm gia eksisoseis loading (0.6*NMCR, 0.96*NMCR)
    
    A = 1.111*rpm^2/n_MCR^2 - 0.0667*rpm/n_MCR;
    
else
    
    A = 1*rpm/n_MCR;

end

B = (1/n_MCR^3*rpm^3);

perithorio = 2-B/A;

end