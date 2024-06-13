% #######################################################################
%  BED TYPE CHARACTERIZATION (Van Rijn 1984)
%  from Bogoni et al. (2017)
% #######################################################################
function   [thetac,fondo]=bed_characterization(Rp,j_thetac,theta,j_bedload)

% sediment size (meters)
      %dsdim = (Rp*Rp / (10^12*1.65*9.81) )^(1/3);

% sediment size (millimeters)
      %dsmm = dsdim*1000;

% dimensionless diameter
      dstar = Rp^(2/3);
      if (dstar<=4)
          thetacVR = 0.24/ dstar;
        elseif ((dstar>4) && (dstar<=10)) 
          thetacVR = 0.14* dstar^(-0.64);
        elseif (dstar>10) && (dstar<=20) 
          thetacVR = 0.04 * dstar^(-0.1);
        elseif ((dstar>20) && (dstar<150))
          thetacVR = 0.013* dstar^0.29;
        elseif (dstar>150)
          thetacVR = 0.055;
      end

% critical Shield number (Brownlie 1981)
      if (j_thetac==1) 
          if j_bedload==0
            thetac = max(0.03,0.5*(0.22*Rp^(-0.6)+0.06*exp(-17.77*Rp^(-0.6)))); %Parker-corrected
          elseif j_bedload==10
            thetac = 0.047; 
          elseif j_bedload==1
            thetac = max(0.0495,0.5*(0.22*Rp^(-0.6)+0.06*exp(-17.77*Rp^(-0.6)))); %Parker-corrected
          elseif j_bedload==2
            thetac = max(0.0386,0.5*(0.22*Rp^(-0.6)+0.06*exp(-17.77*Rp^(-0.6)))); %Parker-corrected
          elseif j_bedload==3
            thetac = max(0.05,0.5*(0.22*Rp^(-0.6)+0.06*exp(-17.77*Rp^(-0.6)))); %Parker-corrected
          end
          
% critical Shield number (Van Rijn 1989)
      elseif (j_thetac==2)  
          thetac=thetacVR;
      end    

% transport parameter
      Tvr = (theta-thetacVR) / thetacVR;

% type classification (Van Rijn 1984)
      if (dstar<=10)

        if (Tvr<=3) 
          fondo = 1;     % plane bed
        elseif (Tvr>3) && (Tvr<=15)
          fondo = 2;     % dune-covered bed
        elseif (Tvr>15)
          fondo = 1;     % plane bed
        end

      elseif (dstar<=100)
        if (Tvr<=0.4)
          fondo=1;       % plane bed
         elseif (Tvr>0.4) && (Tvr<=15)
          fondo = 2;     % dune-covered bed
        elseif (Tvr>15)
          fondo = 1;     % plane bed
        end
      else  
          fondo = 1;     % plane bed
      end 
end
%#######################################################################
%#######################################################################