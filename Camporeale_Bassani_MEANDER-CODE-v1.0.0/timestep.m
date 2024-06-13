function deltat=timestep(Ef,deltas,dU,cstab,tollc,dt0,flag_dt)
%-----------------------------------------------------------------------
% TIME STEP
%-----------------------------------------------------------------------

 year2sec = 86400* 365;

% fixed time step
      if (flag_dt==1) 
        deltat = dt0;
        
% dynamic time step
      elseif (flag_dt==2)

% find maximim value of migration rate
        maxCSI = Ef * max(dU);
      
% find time step threshold (dimensionless years)
        dtmax = cstab * min(deltas,tollc)/(maxCSI) / year2sec;
      
% set current time step
        deltat = min(dt0, dtmax);
        
% wrong flag
      else
        error('ERROR! Wrong flag for time step');
      end
end