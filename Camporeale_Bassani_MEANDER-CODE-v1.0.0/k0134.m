% -----------------------------------------------------------------------
%  Coefficients k0, k1, k3, k4
%  re-adapted from Bogoni et al. (2017)
% -----------------------------------------------------------------------
function [k0,k1,k3,k4]=k0134(Cf,Nz)

g0=zeros(1,Nz);
dg0=g0; g1=g0;dg1=g0; g2=g0;dg2=g0;zg=g0;u0=g0;
GG0=g0; dGG0=g0; GG1=g0; dGG1=g0; 
k0z=g0;k1z=g0;

% -----------------------------------------------------------------------
%  Input data:
      kvk = 0.41; %von karman constant
      A = 1.84;   %coeff. distribution niT
      B = -1.56;  %coeff. distribution niT

      z0 = exp(-kvk/sqrt(Cf) - 0.777);
      deltaz = (1-z0)/(Nz-1);

% -----------------------------------------------------------------------
%  turbulent viscosity nuT0(z) and velocity u0(z)
      I10 = log(z0)-log(1-z0);
      I20 = -(z0+log(1-z0));
      I30 = -(0.5*z0^2 + z0 + log(1-z0));
      const = -(I10 + 2*A*I20 + 3*B*I30);

      xiini = 0;
      xiend = -log(z0);
      deltaxi = (xiend - xiini)/(Nz-1);

% -----------------------------------------------------------------------
%  Computation of G0
% -----------------------------------------------------------------------
      g0(1)  = 0; dg0(1) = 1;g1(1)  = 0;dg1(1) = 1;
      g2(1)  = 0;dg2(1) = 1;

      for j = 1: Nz-1
        xi = xiini + (j-1)*deltaxi;
        zg(j) = z0*exp(xi);
        deltaz = zg(j)*(exp(deltaxi)-1);

% u0(z)
        aux1 = log(zg(j)/z0);
        aux2 = A*(zg(j)^2-z0^2);
        aux3 = B*(zg(j)^3 -z0^3);
        u0(j)= (aux1 + aux2 + aux3)*sqrt(Cf)/kvk;   
   
% g0
        ig=0;zz = zg(j);yj1 = g0(j);
        yj2 = dg0(j);
        [yjj1,yjj2]=Rungeg(A,B,Cf,z0,kvk,0,zz,deltaz,yj1,yj2,ig);
        g0(j+1) = yjj1;  
        dg0(j+1) = yjj2;

% g1
        ig = 1;
        zz = zg(j);
        yj1 = g1(j);
        yj2 = dg1(j);
        [yjj1,yjj2]=Rungeg(A,B,Cf,z0,kvk,0,zz,deltaz,yj1,yj2,ig);
        g1(j+1) = yjj1;  
        dg1(j+1) = yjj2;

% g2
        ig = 2;
        zz = zg(j);
        yj1 = g2(j);
        yj2 = dg2(j);
        [yjj1,yjj2]=Rungeg(A,B,Cf,z0,kvk,0,zz,deltaz,yj1,yj2,ig);
        g2(j+1) = yjj1; 
        dg2(j+1) = yjj2;
      end

      xi = xiini+(Nz-1)*deltaxi;
      zg(Nz) = z0*exp(xi);
      q1 = -dg1(Nz)/dg0(Nz);
      q2 = -dg2(Nz)/dg0(Nz);

% u0(N)
      aux1 = log(zg(Nz)/z0);
      aux2 = A*(zg(Nz)^2 -z0^2);
      aux3 = B*(zg(Nz)^3 -z0^3);
      u0(Nz) = (aux1 + aux2 + aux3)*sqrt(Cf)/kvk;
        
%  g1,g2,g3
     Sumg0=simps(zg,g0);Sumg1=simps(zg,g1);Sumg2=simps(zg,g2);

      ha0 = -(Sumg2 + q2*Sumg0)/(Sumg1 + q1*Sumg0);
      q0 = q1*ha0+q2;

      for j = 1:Nz
        GG0(j) = q0*g0(j) + ha0*g1(j) + g2(j);
        dGG0(j) = q0*dg0(j) + ha0*dg1(j) + dg2(j);
      end

%-----------------------------------------------------------------------
% G1
%-----------------------------------------------------------------------
      g2(1) = 0;      dg2(1) = 1;

      for j = 1: Nz-1
        deltaz = zg(j)*(exp(deltaxi)-1);

% g21
        ig = 21;
        zz = zg(j);
        yj1 = g2(j);
        yj2 = dg2(j);
        GG0aux = (GG0(j)+GG0(j+1))/2;
        [yjj1,yjj2]=Rungeg(A,B,Cf,z0,kvk,GG0aux,zz,deltaz,yj1,yj2,ig);
        g2(j+1) = yjj1;  
        dg2(j+1) = yjj2;
      end

      q2 = -dg2(Nz)/dg0(Nz);

% g1,g2,g3
  
      Sumg0=simps(zg,g0);Sumg1=simps(zg,g1);Sumg2=simps(zg,g2);

      ha1 = -(Sumg2 + q2*Sumg0)/(Sumg1 + q1*Sumg0);
      q0 = q1*ha1+q2;

      for j = 1: Nz
        GG1(j) = q0*g0(j) + ha1*g1(j) + g2(j);
        dGG1(j) = q0*dg0(j) + ha1*dg1(j) + dg2(j);
      end

% -----------------------------------------------------------------------
% k3,k4,k5
% -----------------------------------------------------------------------
      Du0 = sqrt(Cf)*(1/z0 + 2*A*z0 + 3*B*z0^2)/kvk;
      k3 = dGG0(1)/Du0;
      k4 = dGG1(1)/Du0;

% -----------------------------------------------------------------------
% k0,k1,k2
% -----------------------------------------------------------------------
      for j = 1: Nz
        k0z(j) = GG0(j)*u0(j);
        k1z(j) = GG1(j)*u0(j);
      end 
       k0=simps(zg,k0z);k1=simps(zg,k1z);

end
