% ----------------------------------------------------------------------
%  SAVITZKY-GOLAY FILTERING
%  A Savitzky-Golay filter smooth the input signal x of size N into the
%  output signal y, through the following relation:
%  y(i) = 1/h * sum _{j = -(np-1)/2}^{(np-1)/2} aj * x_{i+j}
%  which involves np terms (np must be odd).
% ----------------------------------------------------------------------
function y=savgolfilter_new (np, n, x)

% assign the half-number of points
      nh = round((np-1)/2);
      
% set the Savitzky-Golay coefficients and the normalization term h
      switch (np)
      case(5)
        h = 35;
        a(nh+1:2*nh+1) = [17, 12, -3];
      case(7)
        h = 21;
        a(nh+1:2*nh+1) = [7, 6, 3, -2];
      case(9)
        h = 231;
        a(nh+1:2*nh+1) = [59, 54, 39, 14, -21];
      case(11)
        h = 429;
        a(nh+1:2*nh+1) = [89, 84, 69, 44, 9, -36];
      case(13)
        h = 143;
        a(nh+1:2*nh+1) = [25, 24, 21, 16, 9, 0, -11];
      case(15)
        h = 1105;
        a(nh+1:2*nh+1) = [167, 162, 147, 122, 87,42, -13, -78];
      case(17)
        h = 323;
        a(nh+1:2*nh+1) = [43, 42, 39, 34, 27,18, 7, -6, -21];
      case(19)
        h = 2261;
        a(nh+1:2*nh+1) = [269, 264, 249, 224, 189, 144, 89, 24, -21, -136];
      case(21)
        h = 3059;
        a(nh+1:2*nh+1) = [329, 324, 309, 284, 249, 204,149, 84, 9, -76, -171];
      case(23)
        h = 805;
        a(nh+1:2*nh+1) = [79, 78, 75, 70, 63, 54,43, 30, 15, -2, -21, -42];
      case(25)
        h = 5175;
        a(nh+1:2*nh+1) = [467, 462, 447, 422, 387, 343,...
            287, 222, 147, 62, -33, -138, -253];
      %case default
      %  error('ERROR! Wrong value for Savitzky-Golay filtering');
      end 
      
% set the mirror part of the coefficient array
      for j = 1: nh
        %a(-j) = a(j);
        a(j)=a(j+nh+1);
      end

% set the initial part of the signal equal to the input
      %for i = 1: nh
        y(1: nh) = x(1: nh);
      %end

% filter the central part of the input signal
      for i = nh+1:n-nh
        y(i) = 0;
        for j = -nh: +nh
          y(i) = y(i) + a(j+nh+1) * x(i+j) / h;
        end
      end
      
% set the final part of the signal equal to the input
      %for i = n-nh+1: n
        y(n-nh+1: n) = x(n-nh+1: n);
      %end
      
end
