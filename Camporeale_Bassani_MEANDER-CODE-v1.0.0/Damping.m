function C=Damping(C,jj,ii)

down=zeros(1,ii-jj+1);up=down;
down(1)=C(jj);k=jj+1;kk=2;
while (k<ii)
    down(kk)=(C(k+1)+C(k-1))/2;
    kk=kk+1;
    down(kk)=C(jj+kk-1);
    kk=kk+1;k=kk+jj-1;
end
down(ii-jj+1)=C(ii);

up(1)=C(jj);up(2)=C(jj+1);
kk=3;k=kk+jj-1;
while (k<ii)
    up(kk)=(C(k+1)+C(k-1))/2;
    kk=kk+1;
    up(kk)=C(jj+kk-1);
    kk=kk+1;k=kk+jj-1;
end
up(ii-jj+1)=C(ii);

k=jj;
while (k<ii)
    C(k)=(down(k-jj+1)+up(k-jj+1))/2;  
    k=k+1;
end

end %function