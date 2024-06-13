function [primocut,X,Y,C,s,lim,ncutoff,num]=Cut(i,k,deltas,s,X,Y,C,num,tempo,ncutoff,lc,saveflag,numtot,name1,name2,name3,name4,name5,name6)
% Cutting algorithm

ncutoff=ncutoff+1;
text=['n. cutoff=' num2str(ncutoff)];disp(text);
primocut=1;
lengthcut=abs(s(k)-s(i));
deltaj=k-i-1;numold=num;num=numold-deltaj; 

Xt=[tempo lengthcut lc X(i:i+deltaj)]; %Xt=X of the cut 
Yt=[tempo lengthcut lc Y(i:i+deltaj)];
st=[tempo lengthcut lc s(i:i+deltaj)];
Ct=[tempo lengthcut lc C(i:i+deltaj)];

X(i+1:num)=X(i+1+deltaj:num+deltaj);Y(i+1:num)=Y(i+1+deltaj:num+deltaj);
X(num+1:numold)=NaN;
Y(num+1:numold)=NaN;

if saveflag == 1
    SaveCutoff(Xt,Yt,st,Ct,name1,name2,name3,name4,name5,name6,tempo)
else
end

s=Calcdist1(X,Y,num,numtot);
[C,lim]=Calccurv2(deltas,X,Y,s,num,numtot,primocut);

end 