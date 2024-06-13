function [X,Y,C,s,num,numtot,ind,lim,S,V_mean,v]=Itcurv1(X,Y,C,s,modelflag,b,num,bedload,q,I0,Hmax,dH,d50,deltat,deltas,E,H0,Cf0,numtot,primocut,Delta,Mdat,lung,numeff,saveflag,time,nome8,ncutoff,lapse,savings,mainFolder,qflag)


[C,lim,ind]=Calccurv2(deltas,X,Y,s,num,numtot,primocut);
c=C(1:num)*b;        %dimensionless curvature for HIPS e ZS

% Compute the Sinuosity without the initial straight points
    primoSSS=find(abs(c(1:num))>0.01*mean(abs(c(1:num))),1);
    primoD=find(X(1:num)>H0/(2*Cf0),1);
    primoSS=max(primoSSS,primoD);
    if (isempty(primoSS)==0)
        primoS=primoSS;
        numeff=num-primoS;
    else
        primoS=num-numeff;
    end   
    dist=sqrt((X(num)-X(primoS))^2+(Y(num)-Y(primoS))^2);
    S=(s(num)-s(primoS))/dist;


% --------------------------------------------------------------------------                      
%  update leading parameters
% --------------------------------------------------------------------------
I=I0/S;
[H,~]=Compute_H(q,d50,I,b,Hmax,dH,Delta,bedload); 
ds=d50/H;                           % dimensionless Roughness [-]
beta=b/H;                           % aspect ratio   
theta=I/(Delta*ds);
U=q/(2*b*H);
[trasp,Cf,rpic,F0,CD,CT,phiD,phiT]=resistance(1,bedload,theta,ds,Delta,d50);


s_adim=s(1:num)/b;  %deltas_adim=deltas/b;
if trasp > 0
    if modelflag==1 %HIPS model 
        dU = dU_HIPS(c,num,H,theta,F0,deltas,s_adim,H0,Cf0,Cf,b); 
    elseif modelflag==2 %ZS model 
        [g10,g20,g30,g40,g11,g21,g31,g41,lamb1,lamb2,lamb3,lamb4,Am] = coefZS(beta,theta,rpic,Cf,CD,CT,phiD,phiT,F0,Mdat);
        dU = dU_ZS(c,num,s_adim,g10,g20,g30,g40,g11,g21,g31,g41,lamb1,lamb2,lamb3,lamb4,Am);
        if max(abs(imag(dU)))>1e-13
            error('max(abs(imag(dU)))>1e-13');
        end    
        dU = real(dU);
    end
else 
    dU = zeros(1,num);
end

Ef=E;

% --------------------------------------------------------------------------                      
%  savings
% --------------------------------------------------------------------------
if saveflag == 1 || saveflag == 2
    dU_dim = dU*U;
    v = dU_dim*Ef; %dimensional bank erosion

    if mod(time,lapse) == 0
        %%%%% NOTE the conversion factors to reduce the size of the saved data (from double to int).
        %%%%% During post-processing, they need to be rescaled accordingly.
        X1=(int64(1e4*X(1:num))).';Y1=(int64(1e4*Y(1:num))).';C1=(int64(1e11*C(1:num))).';
        s1=(int64(1e4*s(1:num))).';v1=(int32(1e14*v(1:num))).';
        T = table(X1,Y1,C1,s1,v1);
        nome=['t' num2str(floor(time)) '.csv'];
        cd(savings);writetable(T,nome);cd(mainFolder);
        V_mean=mean(abs(v),'omitnan');
        if saveflag == 1 || saveflag == 2
            pp=load(nome8);
            data.t=[pp.data.t time];
            data.S=[pp.data.S S];data.V_mean=[pp.data.V_mean V_mean];
            Xmin=100;nn=find(X>Xmin); lc=s(num)-s(nn(1));
            data.lc=[pp.data.lc lc];
            if qflag==2, data.q=[pp.data.q q];end
            save(nome8,'-append','data');
        end
        text=['time_save=' num2str(time/365) ' years'];disp(text);
        text=['Cutoffs done=' num2str(ncutoff)];disp(text);
    end
end
% --------------------------------------------------------------------------

% Move the points
[Xn,Yn,V_mean,numtot]=ShiftPoints(X,Y,v,deltat,num,numtot);

X=Xn;Y=Yn;
if X(1)>deltas 
    X=[0 X];Y=[Y(1) Y];C=[0 C];num=num+1;
elseif X(1)<-deltas 
    X=X(2:end); Y=Y(2:end); C=C(2:end);num=num-1;   
end

X=sgolayfilt(X(1:num),3,15); Y=sgolayfilt(Y(1:num),3,15);

s=Calcdist1(X,Y,num,numtot);

[X,Y,s,num,numtot,lim]=SplinePoints(deltas,lim,X,Y,s,num,numtot,lung);
[C,lim,ind]=Calccurv2(deltas,X,Y,s,num,numtot,primocut);


