function [tfin,ro]=StochasticDischarges_generator(Qb,namefile,cv0,tau,dt,base,tempo_end)
% Stochastic discharges generator

%tau = integral scale [days] 
%cv0 = coefficient of variation of Q series

tserie=365*tempo_end; %[days]
cv=(1-2*base/Qb)*cv0; %Coefficient of variation of CPP series

mu=Qb/2-base;         %mean of the Compound Poisson Process (CPP)
alfa=mu*cv^2;
Tm=tau*cv^2; lambda=1/Tm; gammar=1/(mu*cv^2);
eventi=floor(tserie/Tm); 
beta=1/tau;           %reciprocal of the mean interval times among shot events

% Random generations of the jumps and interval times
h=exprnd(alfa,1,eventi); T=exprnd(Tm,1,eventi);dt=0.001;
nt=2; ro1=mu; Ti=0; ro=zeros(1,round(eventi*Tm/dt)); Tf=0;
Delta=zeros(1,eventi);
t=0; ro(1)=mu;
for i=1:eventi
    Tf=Tf+T(i); % exponential decay (deterministic part)
    %tempo=t:dt:Tf;ntt=nt:1:nt+length(tempo)-1;
    %ro(ntt)=ro1*exp(-beta*(tempo-Ti));nt=ntt(end)+1;
     while (t<Tf)
         ro(nt)=ro1*exp(-beta*(t-Ti));
         nt=nt+1;t=t+dt;
     end
    % Jump (shot noise)
    ro1=ro(nt-1)+h(i);
    ro(nt)=ro1;
    t=t+dt;nt=nt+1;Ti=Tf;
end
ro=ro(1:nt-1)+base;  %Stochastic discharges + minimum value
times=dt*(1:nt-1);   %[days]
%f(1)=figure;
%plot(times/365,ro); 
%xlabel('$$years$$','Interpreter','latex'); ylabel('$$Q (m^3/s)$$','Interpreter','latex');
%title('Stochastic discharges')
x=0:0.01:3*mu;
%f(2)=figure;
%plot(x+base,(exp(-x*gammar).*x.^(lambda/beta-1)*gammar^(lambda/beta))/gamma(lambda/beta))
%hold on; histogram(ro,'Normalization','pdf');
%title('Theorethical vs numerical verification')
%
tfin=times(end);
% savefig(f,'TwoFiguresFile.fig')
%
name_Qstoch=['Qstoch_params_' namefile];
save(name_Qstoch,'dt','tau','cv','mu','base','h','T','nt','ro','times','x','gammar','lambda','beta','tfin');

end %function