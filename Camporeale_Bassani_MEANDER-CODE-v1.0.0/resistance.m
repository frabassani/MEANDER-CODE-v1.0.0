% #######################################################################
%  FLOW RESISTANCE AND SEDIMENT TRANSPORT INTENSITY
%  from Bogoni et al. (2017)
% #######################################################################
function [trasp,Cf,rpic,F0,CD,CT,phiD,phiT]=resistance(j_thetac,j_bedload,theta,ds,Delta,d50)


%  j_thetac=1  : thetac              (1=Brownlie 1981; 2=Van Rijn 1989)                            
%  j_bedload   : bed load formula    (see script input_data.m)

flagbed=1;rpic0=0.5;jmodel=1;g=9.81;mi=1e-6; 
% -----------------------------------------------------------------------
%  BED TYPE CHARACTERIZATION
% -----------------------------------------------------------------------
Rp=sqrt(Delta*g*d50^3)/mi;
[thetac,fondo]=bed_characterization(Rp,j_thetac,theta,j_bedload);

% -----------------------------------------------------------------------
%  BED AND TRANSPORT TYPES FROM Rp + VAN RIJN METHOD(1984)
      if (flagbed==1) 

%  FRICTION COEFFICIENT

%  plane bed
        if (fondo==1) 
          [Cf,dCD,dCT,rpic]=resistance_planebed(ds,rpic0);
% dune-covered bed
        elseif (fondo==2)
          [Cf, dCD, dCT,rpic]=resistance_dunebed(ds, theta,  rpic0);
        end

% SEDIMENT TRANSPORT CHARACTERIZATION

% no sediment transport transport
        if (theta<thetac) 
          trasp = 0;
        else
          [trasp,F0]=seditrans_characterization(theta,ds,Rp,Cf);
        end

% SEDIMENT TRANSPORT INTENSITY
        
% bedload
        if (trasp==1)
% Parker 1978
          if (j_bedload==0)
            [phi,dphiD,dphiT]=seditrans_parker78(theta,thetac);
            if phi==0,trasp=0; end

% Meyer-Peter & Muller 1948 modified
          elseif (j_bedload==1)
            [phi,dphiD,dphiT]=seditrans_mpm(theta,thetac);
            if phi==0,trasp=0; end

% Meyer-Peter & Muller 1948 
          elseif (j_bedload==10)
            [phi,dphiD,dphiT]=seditrans_mpm48(theta,thetac);
            if phi==0,trasp=0; end

% Parker (1982, 1990)
          elseif (j_bedload==2)
          [phi, dphiD, dphiT]=seditrans_parker(theta);  
          if phi==0,trasp=0; end

% Ashida & Michiue (1972)        
          elseif (j_bedload==3)
              [phi,dphiD,dphiT]=Ashida_Michiue(theta);
              if phi==0,trasp=0; end
% flag error
          else
            error('ERROR! Wrong flag for bedload transport');
          end

% total load = bedload transport + suspended transport
        elseif (trasp==2)

% Engelund and Hansen (1967)
          [phi, dphiD, dphiT]=seditrans_EH(theta, Cf, dCD, dCT);
          if phi==0,trasp=0; end

        end

% -----------------------------------------------------------------------
%  BED AND TRANSPORT TYPES FROM FROM INPUT SETTING
        elseif (flagbed==2)
%-----------------------------------------------------------------------

% PLANE BED
        if (typebed==1)

% flow resistance and derivatives
          [Cf,dCD,dCT,rpic]=resistance_planebed(ds,rpic0);

% sediment transport intensity and derivatives
% Meyer-Peter & Muller 1948 modified
          if (j_bedload==1)
          [phi,dphiD,dphiT]=seditrans_mpm(theta,thetac);
% Meyer-Peter & Muller 1948
          elseif (j_bedload==10)
          [phi,dphiD,dphiT]=seditrans_mpm48(theta,thetac);
% Parker (1982, 1990)
          elseif (j_bedload==2)
           [phi, dphiD, dphiT]=seditrans_parker(theta);
           if phi==0,trasp=0; end
% flag error
          else
            error('ERROR! Wrong flag for bedload transport');
          end 

% DUNE-COVERED BED
       elseif (typebed==2)

% flow resistance and derivatives
          [Cf, dCD, dCT,rpic]=resistance_dunebed(ds, theta,  rpic0);

% sediment transport intesity and derivatives
% Engelund and Hansen (1967)
          [phi, dphiD, dphiT]=seditrans_EH(theta, Cf, dCD, dCT);
          if phi==0,trasp=0; end

% flag error
        else
         error('ERROR! Wrong flag for bed configuration');
       end
        
% bed characterization
       [trasp,F0]=seditrans_characterization(theta,ds,Rp,Cf);

%-----------------------------------------------------------------------
% end if for resistance and intensity
     end
%-----------------------------------------------------------------------

% parameters for ZS model
      if (jmodel==1 && trasp>0)
        [CD,CT,phiD,phiT]=zs_terms(theta,phi,Cf,dCD,dCT,dphiD,dphiT);
      end
if trasp==0
    CD=0;CT=0; F0=0;phiD=0;phiT=0;
end    

%#######################################################################
% SEDIMENT TRANSPORT TYPE CHARACTERIZATION (Van Rijn 1984)
%#######################################################################
function [trasp,F0]=seditrans_characterization(theta,ds,Rp,Cf)

% sediment size (meters)
      dsdim = (Rp*Rp / (10^12*1.65*9.81) )^(1/3);     

% Froude number
      F0 = sqrt(1.65* ds * theta / Cf);
      F02 = F0^2;

% uniform flow velocity
      U0 = sqrt(F02 * 9.81 * dsdim / ds);

% shear velocity
      ustar = U0 * sqrt(Cf);

% dimensionless settling velocity (Parker 1978, Van Oyen et al. 2011)
      CC = log10(Rp);
      DD = -1.181 + 0.966*CC - 0.1804*CC^2 +0.003746*CC^3 + 0.0008782*CC^4;
      ws0 = 10^DD;

% dimensional settling velocity
      wsdim = ws0 * sqrt(1.65* 9.81* dsdim);

% threshold for suspended load
      Frstar = ustar / wsdim;

% sediment flux characterization (Van Rjin 1984)
      if (Rp>31.62)

% bedload transport
        if (Frstar<0.4)
          trasp = 1;
% bed load + suspended load
        else
          trasp = 2;
        end
      else

% bed load
        if Frstar<(4*Rp^(-2/3))
          trasp = 1;
% bed load + suspended load
        else
          trasp = 2;
        end
      end
end
%#######################################################################
%#######################################################################

% #######################################################################
%  RESISTANCE FOR PLANE BED (Meyer-Peter & Muller, Keulegan 1938)
% #######################################################################
function [Cf,dCD,dCT,rpic]=resistance_planebed(ds,rpic0)

% conductance^-1
      CC = 1/( 6+ 2.5* log(1/(2.5*ds)));

% friction coeffient
      Cf = CC^2;

% derivatives of friction coefficients
      dCD = - 5 * CC^3;
      dCT = 0;
      rpic = rpic0;
end

% #######################################################################
% RESISTANCE FOR DUNE-COVERED BED (Engelund-Hansen 1967)
% #######################################################################
function [Cf, dCD, dCT,rpic]=resistance_dunebed(ds, theta,  rpic0)

% experimental function
      if (theta<0.06)
        theta1 = theta;
        dtheta1 = 1;
      elseif (theta>=0.06) && (theta<=0.55)
        theta1 = 0.06 + 0.4 * theta^2;
        dtheta1 = 0.8 * theta;
%        theta1 = 0.06d0 + 0.3d0 * theta**1.5
%        dtheta1 = 0.45d0 * theta**0.5
      elseif ((theta>0.55) && (theta<=0.8))
        theta1 = 0.06 + 0.4 * theta^2;
        dtheta1 = 0.8 * theta;
%        theta1 = 0.06d0 + 0.3d0 * theta**1.5
%        dtheta1 = 0.45d0 * theta**0.5
      elseif ((theta>0.8) && (theta<1.1068)) 
        A = 0.439218723;
        B = 1.072070082 - A;
        D = 0.85;
        F = 1.1068 - D;
        theta1 = A + B/F * (theta-D);
        dtheta1 = B/F;
      elseif (theta>=1.1068) 
        EE = 0.3 + 0.7*theta^(-1.8);
        theta1 = EE^(-0.56);
        dtheta1 = 0.7056 * theta^(-2.8) * EE^(-1.56);
      end

% Shields ratio
      AA = theta1 / theta;

% derivative of Shields ratio
      BB = dtheta1 / theta - theta1 / theta^2;

% conductance
      CC = 1 / ( 6 + 2.5 * log(AA/(2.5*ds)) );

% friction resistance
      Cf = CC^2 / AA;

% derivatives of friction coefficient
      dCD = - 5 * CC^3 / AA;
      dCT = - CC^2 / AA^2 * BB * (5 * CC + 1);
      rpic = rpic0/sqrt(AA);
%      rpic = rpic0

end
% #######################################################################
%  SEDIMENT TRANSPORT INTENSITY (Parker 1978)
% #######################################################################
function [phi,dphiD,dphiT]=seditrans_parker78(theta,thetacr)

% critical threshold
      %thetacr = 0.03;

% sediment transport intensity
      phi = max(0, 11*theta^(3/2)*(1-thetacr/theta)^(4.5));
     

% derivatives of sediment transport intensity
      dphiD = 0;
      dphiT = max(0, (16.5*(theta - thetacr)^4.5 + 49.5*(theta - thetacr)^3.5*thetacr)/theta^4);
end

% #######################################################################
%  SEDIMENT TRANSPORT INTENSITY (Meyer-Peter & Muller modified)
% #######################################################################
function [phi,dphiD,dphiT]=seditrans_mpm48(theta,thetacr)

% original formula

% critical threshold
%      thetacr = 0.047

% sediment transport intensity
      phi = max(0, 8*(theta-thetacr)^(1.5));

% derivatives of sediment transport intensity
      dphiD = 0;
      dphiT = max(0, 8 * 1.5 *(theta-thetacr)^(1.5-1));
end

% modified formula

% critical threshold
     % thetacr = 0.0495;


% #######################################################################
%  SEDIMENT TRANSPORT INTENSITY (Meyer-Peter & Muller modified)
% #######################################################################
function [phi,dphiD,dphiT]=seditrans_mpm(theta,thetacr)

% original formula

% critical threshold
%      thetacr = 0.047

% sediment transport intensity
%      phi = MAX(0, 8*(theta-thetacr)^(1.5))

% derivatives of sediment transport intensity
%      dphiD = 0
%      dphiT = MAX(0, 8 * 1.5 *(theta-thetacr)^(1.5-1))

% modified formula

% critical threshold
     % thetacr = 0.0495;

% sediment transport intensity
      phi = max(0, 3.97*(theta-thetacr)^(1.5));
     

% derivatives of sediment transport intensity
      dphiD = 0;
      dphiT = max(0, 3.97 * 1.5 *(theta-thetacr)^(1.5-1));
end


% #######################################################################
%  SEDIMENT TRANSPORT INTENSITY (Ashida&Michiue (1972))
% #######################################################################
function [phi,dphiD,dphiT]=Ashida_Michiue(theta)

% critical threshold
      thetacr = 0.05;

% sediment transport intensity
      phi = max(0, 17*(theta-thetacr)*(sqrt(theta)-sqrt(thetacr)));
    

% derivatives of sediment transport intensity
      dphiD = 0;
      dphiT = max(0, 17*(sqrt(theta) - sqrt(thetacr)) + (17*(theta - thetacr))/(2*sqrt(theta)));
end

% #######################################################################
%  SEDIMENT TRANSPORT INTENSITY (Parker 1982, 1990)
% #######################################################################
function [phi, dphiD, dphiT]=seditrans_parker(theta)

      thetar = 0.0386;
      csi = theta / thetar;
      A = 0.00218;

      if (csi>1.65)
        B = 5474;
        F = 0.853;
        D = 4.5;
        G0 = B * (1- F/csi)^D;
        dG = B * D * ( 1- F/csi)^(D-1) * F / (csi^2);
      elseif (csi>=1) && (csi<=1.65)
        B = 14.2;
        F = 9.28;
        G0 = exp(B * (csi-1) - F*(csi-1)^2);
        dG = G0 * (B - 2 * F *(csi-1));
      elseif (csi<1)
        B = 14.2;
        G0 = csi^B;
        dG = B * csi^(B-1);
      end

% sediment transport intensity
      phi = A * G0 * theta^1.5;
   

% derivatives of sediment transport intensity
      dphiD = 0;
      dphiT = A*(dG/thetar)*theta^1.5 + A*G0*1.5*theta^0.5;

end

% #######################################################################
%  SEDIMENT TRANSPORT INTENSITY (Engelund and Hansen 1967)
% #######################################################################
function [phi, dphiD, dphiT]=seditrans_EH(theta, Cf, dCD, dCT)

 A = 0.05;B = 2.50;

% sediment transport intensity
      phi  = A/Cf * theta^B;
      if phi==0,trasp=0; end

% derivatives of sediment transport intensity
      F = -A * theta^B / Cf^2;
      dphiD = F * dCD;
      dphiT = F * dCT + A*B/Cf * theta^(B-1);
end
      
% #######################################################################
%  TERMS FOR ZS MODEL (Zolezzi and Seminara 2001)
% #######################################################################
function [CD,CT,phiD,phiT]=zs_terms(tau,phi,Cf,dCD,dCT,dphD,dphT)

% terms related to the friction coefficient
      CD = dCD / Cf;
      CT = dCT * tau / Cf;
      
% terms related to the sediment intensity
      if (phi>0)
        phiD = dphD / phi;
        phiT = dphT * tau / phi;
      else
        phiD = 0;
        phiT = 0;
      end
end
end