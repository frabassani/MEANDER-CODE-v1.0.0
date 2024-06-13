function [X,Y,s,num,numtot,lim]=SplinePoints(deltas,lim,X,Y,s,num,numtot,lung)

snum=s(num);
if (isnan(snum))                 
    snum_ind = min(find(isnan(s))); 
    snum = s(snum_ind-1);
end
lim=lim-2;

ss(1:num-lim+1)=s(lim:num);  %to neglect a few points at the beginning
x(1:num-lim+1)=X(lim:num);   % "
y(1:num-lim+1)=Y(lim:num);   % "

s(1)=0; 
si=[0:deltas:snum snum];
X1=spline(ss,x,si); Y1=spline(ss,y,si);

num=length(si);
if num>numtot,numtot=numtot+2000;end

X=[X1 NaN(1,numtot-num)]; Y=[Y1 NaN(1,numtot-num)];
if (max(X)>lung)
    [~,i]=find(X>lung);
    X(i(1):end)=NaN;Y(i(1):end)=NaN;num=i(1)-1;
end
s=Calcdist1(X,Y,num,numtot);

end %function