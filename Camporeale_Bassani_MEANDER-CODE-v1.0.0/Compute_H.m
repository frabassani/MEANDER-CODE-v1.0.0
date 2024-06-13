function [h,Cf]=Compute_H(Q,d50,I,b,Hmax,dh,Delta,bedload)
g=9.81;err_rel=zeros(1,1000);Hi=err_rel;Cfi=err_rel;Qi=err_rel;
h=dh; err=0.01;
errore=3*err; j=0;
while errore>err
    h=h+dh; j=j+1;
    ds=d50/h;
    theta=I/(Delta*ds);
    [~,Cf,~,~,~,~,~,~]=resistance(1,bedload,theta,ds,Delta,d50);
    Q1=2*b*h*(Cf^(-1/2))*sqrt(g*h*I);
    err_rel(j)=abs(Q1/Q-1);
    Hi(j)=h;Cfi(j)=Cf;Qi(j)=Q1;
    errore=abs(Q1/Q-1);
    if h>Hmax
        errore=err/2; 
        j=find(err_rel==min(err_rel));
        h=Hi(j);
        Cf=Cfi(j);
    end
end 
    



