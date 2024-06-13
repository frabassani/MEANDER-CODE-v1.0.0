% -----------------------------------------------------------------------
%  ANGLE THETA OF THE RIVER AXIS WITH RESPECT TO THE CARTESIAN REFERENCE
%  from Bogoni et al. (2017)
% -----------------------------------------------------------------------
%  Given a plane curve specified by its cartesian coordinate (x,y), 
%  calculates the angle theta between the local tangent and the x-axis
%  Setting th(j)=-datan2(dy,dx) implies that :
%      c(js) = ( th(js+1)-th(js-1) )/2*ds
%      xnew(js) = x(js) + dt * deltaU(js) * dsin(th(js))
%      ynew(js) = y(js) + dt * deltaU(js) * dcos(th(js))
% -----------------------------------------------------------------------
function th=angle(N,x,y)
th=zeros(1,N);   
% first point
      j = 1;
      dy = y(j+1)-y(j);
      dx = x(j+1)-x(j);
      th(j) = -atan(dy/dx);
      aux0 = 0.5* th(j);
      
% central points      
      for j = 2: N-1
        dy = y(j+1) - y(j);
        dx = x(j+1) - x(j);
        aux1 = -0.5 * atan(dy/dx);
        th(j) = aux1 + aux0;
        aux0 = aux1;
      end 

% last point      
      j = N;
      dy = y(j)-y(j-1);
      dx = x(j)-x(j-1);
      th(j) = -atan(dy/dx);

% set critical angles
      a1 = pi/2;
      a2 = pi*3/2;

% check if th provided sharp gaps crossing pi angles
      j = 1;
      if ((abs(th(j+1)-th(j))>a1) && (abs(th(j+1)-th(j))<a2))     
        th(j) = th(j) - pi * sign(th(j+1)-th(j));
      elseif (abs(th(j+1)-th(j))>pi)
        th(j) = th(j) - 2 * pi * sign(th(j+1)-th(j));
      end

      for j = 2: N
        if ((abs(th(j)-th(j-1))>a1) && (abs(th(j)-th(j-1))<a2))     
          th(j) = th(j) - pi * sign(th(j)-th(j-1));
        elseif (abs(th(j)-th(j-1))> pi) 
          th(j) = th(j) - 2* pi * sign(th(j)-th(j-1));
        end
      end

end
