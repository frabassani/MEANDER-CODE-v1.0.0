function [C,lim,ind]=Calccurv2(deltas,X,Y,s,num,numtot,primocut)
% Computation of initial curvature

ind=0;C=zeros(1,num);
C(1)=0;
for i=2:num-1
    a1=X(i)-X(i-1);a2=Y(i)-Y(i-1);
    b1=X(i+1)-X(i);b2=Y(i+1)-Y(i);
    v1=sqrt((a1^2)+(a2^2));v2=sqrt((b1^2)+(b2^2));
    rapp=((a2*b1)-(a1*b2))/(v1*v2);
    if (rapp>1.0), rapp=1.0; end
    if (rapp<-1.0), rapp=-1.0; end
    C(i)=(asin(rapp))/(s(i)-s(i-1));
    if (abs(C(i))>(1/deltas)), C(i)=sign(C(i))/deltas;end
end
C(num)=0;C(num+1:numtot)=NaN;

if (primocut==1)
    [C,ind]=searchDamping(C,num);
    if (ind==1)
        while (ind==1)
            [C,ind]=searchDamping(C,num);
        end
    end
end
lim=3;
C(num)=0; %last point

end %function