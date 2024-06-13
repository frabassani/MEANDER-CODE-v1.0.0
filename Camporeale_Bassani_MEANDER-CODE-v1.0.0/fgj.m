% --------------------------------------------------------
%  Function for calculating the function fff
% --------------------------------------------------------
function fff=fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,ig)

%--------------------------------------------------------       
      aux1 = log(zrk/z0);
      aux2 = A * (zrk^2 - z0^2);
      aux3 = B * (zrk^3 - z0^3);
      uzero = (aux1 + aux2 +aux3) * sqrt(Cf) / kvk;

      if (zrk<(1- deltaz/4))
        auxnu = 1+ 2*A*(zrk^2) + 3*B*(zrk^3);
        nuTzero = kvk*zrk*(1-zrk)/auxnu;
        Dauxnu = 4*A*zrk+9*B*zrk^2;
        DnuTzero = kvk*((1-2*zrk)*auxnu-zrk*(1-zrk)*Dauxnu)/(auxnu^2);
       else
        nuTzero =   0.0614;
        DnuTzero = -0.0338;
      end

%--------------------------------------------------------       
 switch (ig)

        case(0)
        pf =0;p0 = 0;
        p1 = DnuTzero/nuTzero;

        case(1)
        pf = 1/nuTzero;
        p0 = 0;
        p1 = DnuTzero/nuTzero;

        case(2)
        pf = -(uzero^2)/nuTzero;
        p0 = 0;
        p1 = DnuTzero/nuTzero;

        case(12)
        pf = -1/nuTzero;
        p0 = 0;
        p1 = DnuTzero/nuTzero;

        case(21)
        pf = (uzero*GG0aux)/nuTzero;
        p0 = 0;
        p1 = DnuTzero/nuTzero;

 end
   fff = -p1*yrk2 - p0*yrk1 + pf;   

 end