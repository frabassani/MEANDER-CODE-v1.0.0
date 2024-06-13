function dU=dU_ZS(c,num,s,g10,g20,g30,g40,g11,g21,g31,g41,lamb1,lamb2,lamb3,lamb4,Am)

% -----------------------------------------------------------------------------------------
%  EXCESS NEAR-BANK LONGITUDINAL VELOCITY
%  Zolezzi & Seminara model (2001), with free boundary conditions, i.e., cm1=cm2=cm3=cm4=0.
%  
%  Re-adapted from Bogoni et al. (2017) - parallelized loop
% -----------------------------------------------------------------------------------------


% c      % axis curvature (multiplied by v0)
% num    % number of points (current)

TOLL=1e-04;Mdat=2;

%  UPSTR : upstream propagating influence (i.e. negative-exp convolution)
%  DWSTR : downstream propagating influence (i.e. positive-exp convolution)
%  DWSTR_BC : upstream-propagating donwstream boundary conditions (positive exp)
%  UPSTR_BC : downstream-propagating upstream boundary conditions (negative exp)
%  LOCAL : local effect of the curvature

%  initialize the excess near-bank velocity array
dU = zeros(1,num);
% -----------------------------------------------------------------------
%  loop over Fourier components
for jm = 1: Mdat
    Amjm=Am(jm);
    UPSTR=zeros(1,num);DWSTR=UPSTR;DWSTR_BC=UPSTR;UPSTR_BC=UPSTR;


    %  loop over the eigenvalues
    for jev = 1:4

        % set current eigenvalues
        switch (jev)
            case(1)
                lambda = lamb1(jm);
                gj0 = g10(jm);
                c_bc =0;
            case(2)
                lambda = lamb2(jm);
                gj0 = g20(jm);
                c_bc = 0;
            case(3)
                lambda = lamb3(jm);
                gj0 = g30(jm);
                c_bc = 0;
            case(4)
                lambda = lamb4(jm);
                gj0 = g40(jm);
                c_bc = 0;
        end

        % initialize pointer
        jtoll = 1;

        % UPSTREAM PROPAGATING INFLUENCE (+ DOWNSTREAM BOUNDARY CONDITIONS)
        if (real(lambda)>0)
            for j = 1:num
                if (abs(exp(-(lambda)*(s(j)-s(1))))>TOLL )
                    jtoll = jtoll + 1;
                end
            end



            parfor j = 1:num-1
                jsend = min(j+jtoll-1,num);
                stt=s(j:jsend);ctt=c(j:jsend);
                integrands=exp(-lambda*(stt-s(j))).*ctt;
                CONV=simps(stt,integrands);
                UPSTR(j) = UPSTR(j) - Amjm * gj0 * CONV;
                DWSTR_BC(j) = DWSTR_BC(j)+c_bc*exp(-lambda*(s(num)-s(j)));
            end
            DWSTR_BC(num) = DWSTR_BC(num)+c_bc;

            % DOWNSTREAM PROPAGATING INFLUENCE (+ UPSTREAM BOUNDARY CONDITIONS)
        elseif (real(lambda)<0)
            for j = 1:num
                if (abs(exp(lambda*(s(j)-s(1))))>TOLL)
                    jtoll = jtoll + 1;
                end
            end

            parfor j = 2:num

                jsend = max(j-jtoll+1,1);
                stt=s(jsend:j);ctt=c(jsend:j);
                integrands=exp(lambda*(s(j)-stt)).*ctt;
                CONV=simps(stt,integrands);
                DWSTR(j) = DWSTR(j) + Amjm * gj0 * CONV;

                UPSTR_BC(j) = UPSTR_BC(j) + c_bc * exp(lambda*(s(j)-s(1)));
            end
            UPSTR_BC(1) = UPSTR_BC(1) + c_bc;
        end

        % end loop over the eigenvalues
    end
    % -----------------------------------------------------------------------

    %  local effect of the curvature, longitudinal velocity and excess near-
    %  bank velocity
    g11jm=g11(jm);g21jm=g21(jm);g31jm=g31(jm);g41jm=g41(jm);
    for j = 1:num
        LOCALj = Amjm * c(j) * (g11jm+g21jm+g31jm+g41jm);
        umj = UPSTR(j)+UPSTR_BC(j)+DWSTR(j) + DWSTR_BC(j) + LOCALj;
        dU(j) = dU(j) + 2*umj * (-1)^(jm-1); 
    end
    % end loop over the Fourier components
end
end
