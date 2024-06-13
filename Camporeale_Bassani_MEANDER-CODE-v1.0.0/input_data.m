function nome_conf=input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,bfrac,E,data_name_in,NumSim)
   
% Bedload sediment transport
bedload = 0;    % =0:  Phi=11*tau0^1.5*(1-0.03/tau)^(4.5), Parker (1978)
                % =1:  Phi=3.97*(tau0-0.0495)^(1.5), Meyer-Peter&Muller (1948) modified
                % =10: Phi=8*(tau0-0.047)^(1.5), Meyer-Peter&Muller (1948) 
                % =2:  Phi=0.00218*G0*tau0^(1.5), Parker (1982, 1990)
                % =3:  Phi=17*(tau0-0.05)*(sqrt(tau0)-sqrt(0.05)), Ashida&Michiue (1972)
                
% Variable discharge settings 
nocorr = 1;     % =0: correlated time series (Compound Poissonian process, CPP)
                % =1: white Gamma-distributed noise              
limitflag = 1;  % =0: no bound to discharge
                % =1: upper bound to discharge (i.e., the bankfull value)
                         
% Savings ON/OFF
saveflag = 1;   % =1: ON, save all;
                % =2: ON, save planimetries but not cutoff;
                % =0: OFF, don't save

       
% Simulation end time & saving time
tempo_end = 100000;     %total simulation time in [years] (max for long-tem sim.)
lapse = 10;             %save time in [days]

% Short- or long-term simulations
longshorttermflag = 1;  % =1: short-term sim.
                        % =2: long-term sim.


% -------------------------------------- end of settings --------------------------------------
% ---------------------------------------------------------------------------------------------

if data_name_in == 0 %initial configuration
    longshorttermflag=1;
    cutoffmaxflag=0;
else
    if longshorttermflag == 1

        Ncutoff_max = 30; %maximum number of cutoffs before stopping the simulation

        cutoffmaxflag=1;  %stop after [Ncutoff_max] cutoffs (if also the sinuosity is S > 1.5)
    elseif longshorttermflag == 2
        cutoffmaxflag=2;
    end
end

g = 9.81; rho = 1000;   %gravity [m/s^2]; water density [kg/m^3]
Delta = 1.65;           %relative sediment desity, rho_s/rho [-]
mi=1e-6;                %kinematic viscosity, [m^2/s]  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SIMULATION STARTS
meander;

end