function [Xn,Yn,V_mean,numtot]=ShiftPoints(X,Y,v,deltat_days,num,numtot)

deltat=deltat_days*3600*24; %[seconds]
Xn=X;Yn=Y;
Xn(1)=0;Yn(1)=0;
if num>numtot,numtot=numtot+2000;end
%V_local=NaN(1,num-3);

for k=2:num-2
    b=Y(k+1)-Y(k-1);a=X(k+1)-X(k-1);
    s=sqrt(a^2+b^2);
    Xn(k)=X(k)-v(k)*deltat*(b/s);
    Yn(k)=Y(k)+v(k)*deltat*(a/s);
    %V_local(k)=sqrt((Xn(k)-X(k))^2+(Yn(k)-Y(k))^2);
end
%V_local(V_local>100*mean(V_local))=NaN;

V_mean=mean(abs(v(2:num-2)));
b=2*(Y(2)-Y(1));a=2*(X(2)-X(1));s=sqrt(a^2+b^2);
Xn(1)=X(1)-v(1)*deltat*(b/s);
Yn(1)=Y(1)+v(1)*deltat*(a/s);

b=Y(num)-Y(num-2);a=X(num)-X(num-2);s=sqrt(a^2+b^2);
Xn(num-1)=X(num-1)-v(num-1)*deltat*(b/s);
Yn(num-1)=Y(num-1)+v(num-1)*deltat*(a/s);

% Move the last point
b=Y(num)-Y(num-1);a=X(num)-X(num-1);
s=sqrt(a^2+b^2);
Xn(num)=X(num)-v(num)*deltat*(b/s);
Yn(num)=Y(num)+v(num)*deltat*(a/s);

end %function