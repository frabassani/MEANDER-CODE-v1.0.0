%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     MEANDER CODE
% This is the MAIN SCRIPT TO BE LAUNCHED
% Edited by CARLO CAMPOREALE(1) and FRANCESCA BASSANI(2)
%
% (1) carlo.camporeale@polito.it; (2) francesca.bassani@epfl.ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q0 = 500;      %bankfull discharge defining the river geometry [m^3/s]
d50 = 0.03;    %mean sediment size in [m]
EE = 1e-6;     %Erodibility coefficient [-]
bom = 0.6;     %factor > 0, to reduce the bankfull width [-]; set to 1 otherwise 


modelflag = 1; % =1: HIPS model [Hasegawa (1977) & Ikeda, Parker and Sawai (1981)] 
               % =2: ZS model [Zolezzi and Seminara (2001)]

NumSim = 1;    % set here the NUMBER to be associated with your simulation (e.g., 1)


%%%%%%%%% Simulations with HIPS model
if modelflag == 1
    % --------------------------------------------------------------------------------------------
    % Simulations with constant discharge - initialization
    deltat = 5;                             %timestep of computation in [days]
    disp('HIPS q_const started');
    nsim=['sim' num2str(NumSim) '_start_']; %name of the initial planimetry (sinuosity>1.1)
    cv=0;                                   %coefficient of variation (= 0 for constant discharge)
    
    %%% Initial configuration - the simulation stops at sinuosity S>1.1
    %%% and saves the planimetry (X,Y). To smooth the initial noise.
    %%% Comment the following 2 lines otherwise
    conf_ini = input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,0,NumSim);
    file_in=[conf_ini '.mat']; 
    % --------------------------------------------------------------------------------------------
    
    if exist('file_in','var') ~= 1 %if the initial configuration was not created
        file_in=1;
    else
        file_in=2;
    end

    deltat = 2;                    %timestep of computation in [days]
    
    %%% Simulations HIPS with:
    nsim=['sim' num2str(NumSim) '_']; cv=0; %constant discharge
    input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,file_in,NumSim);
    disp('Starts IPS a');
    nsim=['sim' num2str(NumSim) 'a_']; cv=0.4; %variable discharge with CV=0.4
    input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,file_in,NumSim);
    disp('Starts IPS b');
    nsim=['sim' num2str(NumSim) 'b_']; cv=0.75; %variable discharge with CV=0.75
    input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,file_in,NumSim);



%%%%%%%%% Simulations with ZS model
elseif modelflag == 2
    % --------------------------------------------------------------------------------------------
    % Simulations with constant discharge - initialization
    deltat = 5;                             %timestep of computation in [days]
    disp('ZS q_const started');
    nsim=['sim' num2str(NumSim) '_start_']; %name of the initial planimetry (sinuosity>1.1)
    cv=0;                                   %coefficient of variation (= 0 for constant discharge)
    
    %%% Initial configuration - the simulation stops at sinuosity S>1.1
    %%% and saves the planimetry (X,Y). To smooth the initial noise.
    %%% Comment the following 2 lines otherwise
    conf_ini = input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,0,NumSim);
    file_in=[conf_ini '.mat']; 
    % --------------------------------------------------------------------------------------------
    
    if exist('file_in','var') ~= 1 %if the initial configuration was not created
        file_in=1;
    else
        file_in=2;
    end    

    deltat = 2;                    %timestep of computation in [days]

    %%% Simulations ZS with:
    nsim=['sim' num2str(NumSim) '_']; cv=0; %constant discharge
    input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,file_in,NumSim);
    disp('Starts ZS c');
    nsim=['sim' num2str(NumSim) 'c_']; cv=0.4; %variable discharge with CV=0.4
    input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,file_in,NumSim);
    disp('Starts ZS d');
    nsim=['sim' num2str(NumSim) 'd_']; cv=0.75; %variable discharge with CV=0.75
    input_data(nsim,Q0,d50,modelflag,cv,bom,deltat,0.5,EE,file_in,NumSim);

end %modelflag