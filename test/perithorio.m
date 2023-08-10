function perithorio = perithorio(rpm,n_MCR)
    if rpm>60*2*pi/60 && rpm<95*2*pi/60
        A = 1.111*rpm^2/n_MCR^2-0.0667*rpm/n_MCR;
    else
        A = 1*rpm/n_MCR;
    end
  %  B = (1/n_MCR^3*rpm^3);
      B = (1/n_MCR^3*rpm^3);

    %perithorio = ((1.111*rpm^2/n_MCR^2-0.067*rpm/n_MCR)-1/n_MCR^3*rpm^3)/(1/n_MCR^3*rpm^3)+1;
  %  perithorio = 1-A/B;
        perithorio = 2-B/A;


end