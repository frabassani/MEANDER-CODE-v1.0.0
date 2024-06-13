% ---------------------------------------------------------------------
%  Runge-Kutta for second order ordinary differential equation
%  input: zz,yj1,yj2       output yjj1,yjj2
% ---------------------------------------------------------------------
 function [yjj1,yjj2]=Rungeg(A,B,Cf,z0,kvk,GG0aux,zz,deltaz,yj1,yj2,ig)

%---------------------------------------------------------------------
      zrk=zz;      yrk1=yj1;      yrk2=yj2;
      fff=fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,ig);
      k11=deltaz*yrk2;
      k12=deltaz*fff;
  
      zrk=zz+deltaz/2;
      yrk1=yj1+k11/2;
      yrk2=yj2+k12/2;
      fff=fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,ig);
      k21=deltaz*yrk2;
      k22=deltaz*fff;
   
      zrk=zz+deltaz/2;
      yrk1=yj1+k21/2;
      yrk2=yj2+k22/2;
      fff=fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,ig);
      k31=deltaz*yrk2;
      k32=deltaz*fff;
   
      zrk=zz+deltaz;
      yrk1=yj1+k31;
      yrk2=yj2+k32;
      fff=fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,ig);
      k41=deltaz*yrk2;
      k42=deltaz*fff;

      yjj1=yj1+(k11 + 2*(k21+k31)+k41)/6;
      yjj2=yj2+(k12 + 2*(k22+k32)+k42)/6;

%---------------------------------------------------------------------
end

