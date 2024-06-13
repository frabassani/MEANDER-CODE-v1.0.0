function [g10,g20,g30,g40,g11,g21,g31,g41,lamb1,lamb2,lamb3,lamb4,Am]=coefZS(beta,theta,rpic,Cf,CD,CT,phiD,phiT,F0,Mdat)

%-----------------------------------------------------------------------
% Coefficients for the Flow Field with Variable Curvature 
% under Periodic Boundary Conditions (Zolezzi and Seminara 2001)
%-----------------------------------------------------------------------

theta0=theta;
coef=zeros(1,5);wj=zeros(1,5);
% -----------------------------------------------------------------------
%  coefficients k0, k1, k2, k3, k4
% -----------------------------------------------------------------------
[k0,k1,k3,k4]=k0134(Cf,1000);
% -----------------------------------------------------------------------
%  ZS coefficients
% -----------------------------------------------------------------------

      s1 = 2 / (1-CT);
      s2 = CD / (1-CT);
      f1 = 2 * phiT / (1-CT);
      f2 = phiD + CD * phiT / (1-CT);

% coefficients a1 - a6
      a1 = beta * Cf * s1;
      a2 = beta * Cf * (s2 - 1);
      a3 = beta * Cf;
      a4 = f1;
      a5 = f2;
      a6 = rpic/(beta*theta0^0.5);

% coefficients b1 - b6
      b1 = - beta*Cf;
      b2 = 1 - sqrt(Cf) * k3;
      b3 = -k0 / (beta*Cf^0.5) - k4/beta;
      b4 = k3 * theta0^0.5 / (rpic*Cf^0.5);
      b5 = - k1 / (Cf*beta^2);
      b6 = k4 * theta0^0.5 / (beta*Cf*rpic);

% coefficients hbar and dbar
      h1bar = b2;
      h2bar = b3;
      h3bar = b5;
      d1bar = (F0^2) * h1bar - b4;
      d2bar = (F0^2) * h2bar - b6;
      d3bar = (F0^2) * h3bar;

      alf0 = a2;
      alf1 = 1 / (F0^2);

      bet2 = a1;
      bet3 = 1;

      gam2 = b1 - a2 * d1bar;
      gam3 = -h1bar - a2 * d2bar;

      del1 = a5 - 1 - (F0^2) * a3 * a6;
      del2 = - (F0^2) * a6;

      eps3 = a4 - 1 - (F0^2) * a3 * a6;
      eps4 = del2;

      et3 = - del1 * d1bar;
      et4 = - del1 * d2bar + (F0^2.d0) * a6 * d1bar;
      et5 = - del1 * d3bar + (F0^2.d0) * a6 * d2bar;
      %et6 = (F0^2) * a6 * d3bar;

% -----------------------------------------------------------------------
Am=zeros(1,Mdat);lamb1=Am;lamb2=Am;lamb3=Am;lamb4=Am;
g10=Am;g20=Am;g30=Am;g40=Am;
g11=Am;g21=Am;g31=Am;g41=Am;
     for jm = 1: Mdat
        M = ( 2*(jm-1) + 1) * pi / 2;
        Am(jm) = ((-1)^(jm-1)) * 2/(M^2);
   
        alf2 = (1-a5) / ( (M^2.d0)*(F0^2)*a6);
 
        bet4  =(1-a4) / ( (M^2.d0)*(F0^2)*a6);
   
        gam4 = - alf2 * d1bar - h2bar - a2 * d3bar;
        gam5 = - alf2 * d2bar - h3bar;
   
        del0 = - (M^2) * a6;
   
        Deltaa = del2 * alf1 - del1 * alf2;
        Delta0 = (del2 * alf0 - del0 * alf2) / Deltaa;
        Delta1 = del2 * Delta0 - del1;
        Delta2 = Delta1 * Delta0 + del0;

        T1 = - del2 * bet2 / Deltaa;
        T2 = - (del2 * bet3 - alf2 * eps3) / Deltaa;
        T3 = - (del2 * bet4 - alf2 * eps4) / Deltaa;

        Tc1 = del2 * gam2 / Deltaa;
        Tc2 = (del2 * gam3 - alf2 * et3) / Deltaa;
        Tc3 = (del2 * gam4 - alf2 * et4) / Deltaa;
        Tc4 = (del2 * gam5 - alf2 * et5) / Deltaa;
%       Tc5=0
         
        csi1 = -Deltaa * Delta1 * T1;
        csi2 = Deltaa * (-Delta1 * T2 + del2 * T1 + eps3);
        csi3 = Deltaa * (-Delta1 * T3 + del2 * T2 + eps4);
        csi4 = Deltaa * del2 * T3;

        mu1 = Deltaa * Delta1 * Tc1;
        mu2 = Deltaa * (Delta1 * Tc2 - del2 * Tc1 + et3);
        mu3 = Deltaa * (Delta1 * Tc3 - del2 * Tc2 + et4);
        mu4 = Deltaa * (Delta1 * Tc4 - del2 * Tc3 + et5) ;
%       mu5=0   
%       mu6=0

        sigma0 = (Delta0 * csi1 + Deltaa * Delta2 * T1) / csi4;
        sigma1 = (csi1 + Delta0 * csi2 + Deltaa * Delta2 * T2) / csi4;
        sigma2 = (csi2 + Delta0 * csi3 + Deltaa * Delta2 * T3) / csi4;
        sigma3 = (csi3 + Delta0 * csi4) / csi4;

        rho0 = (Delta0 * mu1 - Deltaa * Delta2 * Tc1) / csi4;
        rho1 = (mu1 + Delta0 * mu2 - Deltaa * Delta2 * Tc2) / csi4;
        rho2 = (mu2 + Delta0 * mu3 - Deltaa * Delta2 * Tc3) / csi4;
        rho3 = (mu3 + Delta0 * mu4 - Deltaa * Delta2 * Tc4) / csi4;
        rho4 =  mu4 /csi4;
%       rho5=0
%       rho6=0
% 
% -----------------------------------------------------------------------
%  quartic roots
% -----------------------------------------------------------------------
        coef(1) = sigma0;
        coef(2) = sigma1;
        coef(3) = sigma2;
        coef(4) = sigma3;
        coef(5) = 1;     

        rooots=roots(fliplr(coef));
        
        for i = 1:3
            for j = i+1:4
               if (rooots(j)>rooots(i)) 
                   x = rooots(i);rooots(i) = rooots(j);
                   rooots(j) = x;
               end 
            end 
        end 
        %%%%%%%%%%%

        lamb1(jm) = rooots(1);
        lamb2(jm) = rooots(2);
        lamb3(jm) = rooots(3);
        lamb4(jm) = rooots(4);

% -----------------------------------------------------------------------
%  Wronskians
% -----------------------------------------------------------------------
        lam1 = rooots(1);
        lam2 = rooots(2);
        lam3 = rooots(3);
        lam4 = rooots(4);

        waux1 = lam3 * lam4 * (lam4-lam3);
        waux2 =-lam2 * lam4 * (lam4-lam2);
        waux3 = lam2 * lam3 * (lam3-lam2);
        W1det =-(waux1 + waux2 + waux3);

        wdet = -lam2 * lam3 * lam4 * W1det;

        waux1 = lam3 * lam4 * (lam4-lam3);
        waux2 =-lam1 * lam4 * (lam4-lam1);
        waux3 = lam1 * lam3 * (lam3-lam1);
        W2det = (waux1 + waux2 + waux3);
                
        wdet = wdet - lam1 * lam3 * lam4 * W2det;

        waux1 = lam2 * lam4 * (lam4-lam2);
        waux2 =-lam1 * lam4 * (lam4-lam1);
        waux3 = lam1 * lam2 * (lam2-lam1);
        W3det =-(waux1 + waux2 + waux3);

        wdet = wdet - lam1 * lam2 * lam4 * W3det;

        waux1 = lam2 * lam3 * (lam3-lam2);
        waux2 =-lam1 * lam3 * (lam3-lam1);
        waux3 = lam1 * lam2 * (lam2-lam1);
        W4det = (waux1 + waux2 + waux3);

        wdet = wdet - lam1 * lam2 * lam3 * W4det;

        wj(1) = W1det / wdet;
        wj(2) = W2det / wdet;
        wj(3) = W3det / wdet;
        wj(4) = W4det / wdet;
% 
% -----------------------------------------------------------------------
%  gj0
% -----------------------------------------------------------------------
        g10(jm) = wj(1) * ( rho0 + rho1*lam1...
              + rho2*lam1^2 + rho3*lam1^3 + rho4*lam1^4 );

        g20(jm) = wj(2) * ( rho0 + rho1*lam2 ...
             + rho2*lam2^2 + rho3*lam2^3 + rho4*lam2^4);

        g30(jm) = wj(3) * ( rho0 + rho1*lam3 ...
               + rho2*lam3^2 + rho3*lam3^3 + rho4*lam3^4);

        g40(jm) = wj(4) * ( rho0 + rho1*lam4 ...
               + rho2*lam4^2+ rho3*lam4^3 + rho4*lam4^4);

% -----------------------------------------------------------------------
%  gj1
% -----------------------------------------------------------------------
        g11(jm) = wj(1) * ( rho1 + rho2*lam1 + rho3*lam1^2 + rho4*lam1^3 );

        g21(jm) = wj(2) * ( rho1 + rho2*lam2 + rho3*lam2^2 + rho4*lam2^3);

        g31(jm) = wj(3) * ( rho1 + rho2*lam3 + rho3*lam3^2 + rho4*lam3^3 );

        g41(jm) = wj(4) * ( rho1 + rho2*lam4 + rho3*lam4^2 + rho4*lam4^3);

      end
