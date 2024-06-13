function [C,Ind]=searchDamping(C,num)

i=2;Ind=0;
while (i<num-10)&&(Ind~=1)
    d0=C(i);d1=C(i)-C(i-1);
    d2=C(i+1)-C(i);d3=C(i+2)-C(i+1);

    if (((d1*d2)<0)&&((d3*d2)<0))
        j=i;i=i+1;n=1;
        d1=C(i)-C(i-1);d2=C(i+1)-C(i);d3=C(i+2)-C(i+1);d4=C(i+3)-C(i+1);
        while (i<num-10)&&((((d1*d2)<0.0)&&((d3*d2)<0.0))||(((d1*d2)<0.0)&&((d2*d3)>0.0)&&((d4*d3)<0.0))||(((d1*d2)>0.0)&&((d2*d3)<0.0)&&((d4*d3)<0.0)))
            n=n+1;i=i+1;
            d1=C(i)-C(i-1);d2=C(i+1)-C(i);d3=C(i+2)-C(i+1);d4=C(i+3)-C(i+2);
        end
        if 	(n>5)
            ii=i+5;jj=j-5;
            if jj<1,jj=1;end
            C=Damping(C,jj,ii);Ind=1;
        end
    end
    i=i+1;
end

end