function perithorio = perithorio(rpm,n_MCR)
    A = 1.111*rpm^2/n_MCR^2-0.067*rpm/n_MCR;
  %  B = (1/n_MCR^3*rpm^3);
      B = (1/n_MCR^3*rpm^3);

    %perithorio = ((1.111*rpm^2/n_MCR^2-0.067*rpm/n_MCR)-1/n_MCR^3*rpm^3)/(1/n_MCR^3*rpm^3)+1;
  %  perithorio = 1-A/B;
        perithorio = 2-B/A;


end