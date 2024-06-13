function s=Calcdist1(X,Y,num,numtot)
%Computation of curvilinear coordinate s

deltas_vec=zeros(1,num);
s=zeros(1,numtot); 
for i=2:num
	deltas_vec(i)=sqrt(((X(i)-X(i-1))^2)+((Y(i)-Y(i-1))^2));
	s(i)=s(i-1)+deltas_vec(i);
end
s(num+1:numtot)=NaN;

end %function
