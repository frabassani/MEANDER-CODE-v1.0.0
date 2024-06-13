% -----------------------------------------------------------------------
%  DISCRETE DERIVATIVE OF A FUNCTION (x,y) OF N POINTS
%  from Bogoni et al. (2017)
% -----------------------------------------------------------------------
%  Flag is a switch for the numerical differentiation method:
%  1 = derivative in a point is computed by averaging the slopes of the 
%      two segments sharing that points
%  2 = derivative is computed with second-order central finite difference
%  3 = derivative is computed with quasi fourth-order central finite difference
%  4 = derivative is computed with fourth-order central finite difference
%      (last solutions requires N >=6 )
%  NB requires the function slope (reported below), to check eventual 
%  denominator approaching zero
% -----------------------------------------------------------------------
function dydx=derivative(N, x, y, flag)

dydx=zeros(1,N);

%-----------------------------------------------------------------------
% BACK-AND-FORTH SLOPE WEIGHTING
      if flag==1 
%-----------------------------------------------------------------------

% forward finite difference for first point
        i = 1;
        dy = y(i+1) - y(i);
        dx = x(i+1) - x(i);
        aux1 = slope(dx,dy);
        dydx(i) = aux1;

% back and forth slope averaging for central points
        for i = 2: N-1
          dy = y(i+1) - y(i);
          dx = x(i+1) - x(i);
          aux2 = slope(dx,dy);

          dydx(i) = 0.5 * (aux1 + aux2);
          aux1 = aux2;
        end 

% backward finite difference for first point
        i = N;
        dy = y(i) - y(i-1);
        dx = x(i) - x(i-1);
        dydx(i) = slope(dx,dy);

%-----------------------------------------------------------------------
% SECOND-ORDER CENTRAL FINITE DIFFERENCE
      elseif flag==2 
%-----------------------------------------------------------------------

% forward finite difference for first point
        i = 1;
        dy = y(i+1) - y(i);
        dx = x(i+1) - x(i);
        dydx(i) = slope(dx,dy);

% central finite difference for central points
        for i = 2: N-1
          dy = y(i+1)-y(i-1);
          dx = x(i+1)-x(i-1);
          dydx(i) = slope(dx,dy);
        end 

% backward finite difference for last point
        i = N;
        dy = y(i) - y(i-1);
        dx = x(i) - x(i-1);
        dydx(i) = slope(dx,dy);

%-----------------------------------------------------------------------
% QUASI FOURTH-ORDER FINITE DIFFERENCE
      elseif flag==3 
%-----------------------------------------------------------------------

% first point
        i = 1;
        dy = y(i+1) - y(i);
        dx = x(i+1) - x(i);
        dydx(i) = slope(dx,dy);

% second point
        i = 2;
        dy = y(i+1) - y(i-1);
        dx = x(i+1) - x(i-1);
        dydx(i) = slope(dx,dy);
      
% central finite difference for central points
        for i = 3: N-2
          dy = -1/12*y(i+2) + 2/3*y(i+1)-2/3*y(i-1) + 1/12*y(i-2);
          dx = (x(i+1)-x(i-1))/2;
          dydx(i) = slope(dx,dy);
        end

% second to last point
        i = N-1;
        dy = y(i+1) - y(i-1);
        dx = x(i+1) - x(i-1);
        dydx(i) = slope(dx,dy);

% last point
        i = N;
        dy = y(i) - y(i-1);
        dx = x(i) - x(i-1);
        dydx(i) = slope(dx,dy);

%-----------------------------------------------------------------------
% FOURTH-ORDER FINITE DIFFERENCE
      elseif (flag==4) 
%-----------------------------------------------------------------------

% forward finite difference for first and second points
        for i = 1: 2
          dy = - 3/2*y(i+2) + 2*y(i+1) - 1/2*y(i);
          dx = (x(i+1)-x(i));
          dydx(i) = slope(dx,dy);
        end 

% central finite difference for central points
        for i = 3: N-2
          dy = -1/12*y(i+2) + 2/3*y(i+1)-2/3*y(i-1) + 1/12*y(i-2);
          dx = (x(i+1)-x(i-1))/2;
          dydx(i) = slope(dx,dy);
        end

% backward finite difference for second-last and last points
        for i = N-1:N
          dy = 1/2*y(i-2) - 2*y(i-1) + 3/2*y(i);
          dx = (x(i)-x(i-1))/2;
          dydx(i) = slope(dx,dy);
        end

      else
        message('ERROR! Wrong flag for derivative subroutine');
      end
end

%-----------------------------------------------------------------------
% SLOPE CHECKING FOR EVENTUAL ZERO DENOMINATOR
%-----------------------------------------------------------------------
function sl=slope(dx,dy)

toll=1e-6;

% if the denominator approaches zero, then set the ratio equal to the
% largest available finite number
      if (abs(dx)<toll)
        sl = realmax * sign(dy) * sign(dx);
      else
        sl = dy/dx;
      end
end
