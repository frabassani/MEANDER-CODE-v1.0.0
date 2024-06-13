function dU=dU_HIPS(C,num,H,theta,F0,deltas,s,H0,Cf0,Cf,b) 

% -----------------------------------------------------------------------------------------
%  EXCESS NEAR-BANK LONGITUDINAL VELOCITY
%  Hasegawa (1977) - Ikeda, Parker & Sawai (1981) model
%  
%  Re-adapted from Camporeale et al. (2005) - parallelized loop
% -----------------------------------------------------------------------------------------


chi1=0.077/(sqrt(Cf)); chi=chi1-(1/3);
beta_=0.55/sqrt(theta);
A=(1./(0.005929*beta_))*(2/45)*((chi+(2/7))/(chi+(1/3)))-1; %slope factor
P=(F0^2+A)/2; % Parker's number
p0=2*Cf*b/H;

span=s(1:num);
C_active=C(1:num);
dU = zeros(1,num);

ismax1=round(50*H0/(2*Cf0)/deltas);
parfor is = 2:ismax1
    indt=1:is; t=span(indt);
    integrands=C_active(indt).*exp(-p0*(span(is)-t));
    dU(is)=2*(-C_active(is)+p0*(P+1)*trapz(t,integrands));
end
parfor is = ismax1+1:length(span(1:end))
    indt=1:is; indt=indt(end-ismax1+1:end); t=span(indt);
    integrands=C_active(indt).*exp(-p0*(span(is)-t));
    dU(is)=2*(-C_active(is)+p0*(P+1)*trapz(t,integrands));
end

end