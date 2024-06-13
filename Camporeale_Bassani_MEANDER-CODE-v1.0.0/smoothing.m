% -----------------------------------------------------------------------
%  CURVATURE SMOOTHING
% -----------------------------------------------------------------------
function c_out=smoothing(N, c)

      aux= c;

      for j = 2: N-1
	    aux(j) = ( c(j-1) + 2*c(j) + c(j+1) ) / 4;
      end 

      c_out= aux;

end