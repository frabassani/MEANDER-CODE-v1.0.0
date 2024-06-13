function [X,Y,C,s,lim,IND0,num,primocut,ncutoff]=Control_Cutoff(C,b,estremi,numcell,numcol,numrig,ptimax,deltas,dim,X,Y,s,num,primocut,tempo,lim,ncutoff,lc,saveflag,numtot,nome1,nome2,nome3,nome4,nome5,nome6)

% Algorithm for filling the address matrix A(numcell,ptimax) 
% Note: The emergency matrix EM can contain a maximum of 10 oversampled cut-offs.
% The cell width is equal to the river width divided by the square root of two.
% from Camporeale et al. (2005)


x=X(1:num);y=Y(1:num);
A=zeros(numcell,ptimax);N=zeros(1,numcell);EM=zeros(10,ptimax+1);
RC=zeros(numcol,2);
RC(:,1)=numrig;RC(:,2)=1;
ylim=(((numrig*dim)-(estremi(4)-estremi(3)))/2)+estremi(4);

for i=2:num
    c=(floor((x(i)-estremi(1))/dim))+1;
    r=(floor((ylim-y(i))/dim))+1;
    rr=(numrig*(c-1))+r;
    if (RC(c,1)>r), RC(c,1)=r;end
    if (RC(c,2)<r), RC(c,2)=r;end
    if (N(rr)<ptimax)
        N(rr)=N(rr)+1;pos=N(rr);A(rr,pos)=i;
    else
        N(rr)=N(rr)+1;ind=0;ind1=0;ind2=0;
        for h=1:10
            if (EM(h,1)~=0), ind=ind+1;end
            if (EM(h,1)==rr)
                hh=2;
                while ((hh<=ptimax+1)&&(ind2==0))
                    if (EM(h,hh)==0),EM(h,hh)=i;ind2=1;end
                    hh=hh+1;
                end
                ind1=1;
            end
        end
        if (ind1==1),EM(ind+1,1)=rr;EM(ind+1,2)=i;end
    end
end


% Cutoff search
IND0=0;c=1;
while ((c<=numcol)&&(IND0~=1))
    r1=RC(c,1);r2=RC(c,2);r=r1;
    while ((r<=r2)&&(IND0~=1))
        rr=(numrig*(c-1))+r;
        if(N(rr)>ptimax)
            h=1;
            while ((h<=10)&&(IND0~=1))
                if (EM(h,1)==rr)
                    pos2=(N(rr)-ptimax)+1;i=A(rr,1);k=EM(h,pos2);
          
                    [primocut,X,Y,C,s,lim,ncutoff,num]=Cut(i,k,deltas,s,X,Y,C,num,tempo,ncutoff,lc,saveflag,numtot,nome1,nome2,nome3,nome4,nome5,nome6);
                    x=X(1:num);y=Y(1:num);
                    IND0=1;
                end
                h=h+1;
            end
        end
        if ((N(rr)<=ptimax)&&(N(rr)~=0)&&(IND0~=1))
            pos=N(rr);
            if (A(rr,pos)>=(A(rr,1)+(3*ptimax)))
                i=A(rr,1);k=A(rr,pos);
                [primocut,X,Y,C,s,lim,ncutoff,num]=Cut(i,k,deltas,s,X,Y,C,num,tempo,ncutoff,lc,saveflag,numtot,nome1,nome2,nome3,nome4,nome5,nome6);
                x=X(1:num);y=Y(1:num);
                IND0=1;
            else
                an=zeros(1,4);
                if ((((rr+numrig-1)<=numcell)&&(N(rr+numrig-1)~=0))&&((abs(A(rr+numrig-1,1)-A(rr,pos)))>(3*ptimax))), an(1)=rr+numrig-1;end
                if ((((rr+numrig)<=numcell)&&(N(rr+numrig)~=0))&&((abs(A(rr+numrig,1)-A(rr,pos)))>(3*ptimax))), an(2)=rr+numrig;end
                if ((((rr+numrig+1)<=numcell)&&(N(rr+numrig+1))~=0)&&((abs(A(rr+numrig+1,1)-A(rr,pos)))>(3*ptimax))), an(3)=rr+numrig+1;end
                if ((((rr+1)<=numcell)&&(N(rr+1)~=0))&&((abs(A(rr+1,1)-A(rr,pos)))>(3*ptimax))), an(4)=rr+1;end
                nn=1;
                while ((nn<=4)&&(IND0~=1))
                    rif=an(nn);
                    if (rif~=0), pos2=N(rif);end
                    if ((rif~=0)&&(pos2>0)&&(pos2<=ptimax))
                        ii=1;
                        while ((ii<=pos)&&(IND0~=1))
                            jj=pos2;
                            while ((jj>=1)&&(IND0~=1))
                                i=A(rr,ii);j=A(rif,jj);
                                dist=sqrt(((x(i)-x(j))^2)+((y(i)-y(j))^2));
                                if (dist<=b)&&((abs(i-j)>(3*ptimax)))
                                    if (i<j)                                      
                                        [primocut,X,Y,C,s,lim,ncutoff,num]=Cut(i,j,deltas,s,X,Y,C,num,tempo,ncutoff,lc,saveflag,numtot,nome1,nome2,nome3,nome4,nome5,nome6);
                                        x=X(1:num);y=Y(1:num);
                                    else                                      
                                        [primocut,X,Y,C,s,lim,ncutoff,num]=Cut(j,i,deltas,s,X,Y,C,num,tempo,ncutoff,lc,saveflag,numtot,nome1,nome2,nome3,nome4,nome5,nome6);
                                        x=X(1:num);y=Y(1:num);
                                    end
                                    IND0=1;
                                end
                                jj=jj-1;
                            end
                            ii=ii+1;
                        end
                    end
                    nn=nn+1;
                end
            end
        end
        r=r+1;
    end
    c=c+1;
end
end 