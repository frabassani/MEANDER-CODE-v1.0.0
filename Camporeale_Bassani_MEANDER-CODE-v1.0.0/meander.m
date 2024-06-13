%%%%%%%%% SIMULATION STARTS

% Discharge (auto-)setting: =1: constant flow;  =2: variable flows
if cv>0 %Q=variable
    qflag = 2; 
    qStr = ['v_cv' num2str(100*cv)];
else    %Q=constant
    qflag = 1; 
    qStr = 'q'; 
end
qBaseStr = 'q';

if modelflag == 1 %HIPS model
    mStr = 'HIPS';
else              %ZS model
    mStr = 'ZS';
end

namefile=[nsim qStr '_' mStr];
mkdir(namefile); %creates folder for saving
mainFolder=pwd;
savings=[mainFolder filesep namefile];

t=tempo_end*365; %[days]

% Hydraulic geometry formulas (Parker et al. (2007)):
Bbf=bom*4.63*(Q0^0.4)/(g^(1/5))*(Q0/sqrt(g*d50^5))^0.0667; %bankfull width [m^3/s] times the reducing factor (set in the main file)
Hbf=0.382*(Q0^0.4)/(g^(1/5));                              %bankfull depth [m]
Jbf=0.101*(Q0/sqrt(g*d50^5))^(-0.344);                     %bankfull slope 

% Parameters computed form INPUT data:
I0=Jbf*5;                          %slope of straight channel [-]
b=0.5*Bbf;                         %half channel width [m]
deltas=b*bfrac;                    %fitting step [m]
dim=3*b/(sqrt(2));
H0=Hbf;
Hmax=H0*10; dH=H0/1000;
[H0,~]=Compute_H(Q0,d50,I0,b,...
    Hmax,dH,Delta,bedload);        %water depth [m] and friction factor Cf (uniform flow conditions), first trial
U0=Q0/(2*b*H0);                    %stream velocity (uniform flow conditions) [m/s]
ds=d50/H0;                         %dimensionless Roughness [-]
beta=b/H0;                         %aspect ratio
tau0=rho*g*H0*I0;                  %bed shear stress of normal flow, [N/m^2]
theta=I0/(Delta*ds);               
[~,Cf0,~,~,~,~,~,~]=resistance(1,bedload,theta,ds,Delta,d50);
Rp=sqrt(Delta*g*d50^3)/mi;
[thetac,~]=bed_characterization(Rp,1,theta,bedload);
Hc=thetac*d50*Delta/Jbf;
Qtrial=Q0/1000;
[Htrial,~]=Compute_H(Qtrial,d50,Jbf,b,Hmax,dH,Delta,bedload);
while Qtrial<Q0 && Htrial<Hc
    [Htrial,~]=Compute_H(Qtrial,d50,Jbf,b,Hmax,dH,Delta,bedload);
    Qtrial=Qtrial+Q0/1000;
end
base=Qtrial;


D0=H0/(2*Cf0); T0=(D0^2/(b*E*U0))/(3600*24*365);
lung=floor(50*13*D0);                             %from Camporeale et al. (2005)
numtot=lung+1000;


% Setting for cutoff search
ncutoff=0; l=0; scartato=1; lim=3; IND0=1;
primocut=1; % =0 if no cutoff yet

% Initial configuration
if data_name_in==0 || data_name_in==1 
    X=0:deltas:lung;
    Y=(b/32)*rand(1,lung)-b/64; Y(1:200)=0;
    nome_conf=['conf_in_' namefile]; save(nome_conf,'X','Y');
    time0=0;Y=[zeros(1,200) Y];
    num=length(X); numeff=num;
else
    data_name_in=['conf_in_sim' num2str(NumSim) '_start_' qBaseStr '_' mStr '.mat'];
    data_ing=load(data_name_in); X=(data_ing.X).'; Y=(data_ing.Y).';
    clear  data;
    time0=0;
    pp=find(isnan(X)==1); 
    if isempty(pp) 
        num=length(X);
    else
        num=min(pp)-1;
    end
    numtot=num+2000;numeff=num;   
end


if (qflag == 2)  %Q var
    if nocorr == 0
        [t,Qstoch]=StochasticDischarges_generator(Q0,namefile,cv,60,deltat,base,tempo_end);
    else
        mu=Q0-base;
        aa=1/cv^2; bb=mu*cv^2;
        pd=makedist('Gamma','a',aa,'b',bb);
        sz=floor(t/deltat);
        Qstoch=random(pd,[1,sz])+base;
    end
end


numt=floor(t/deltat); Xmin=100;

text=['Simulation starts ' namefile];
disp(text);

% Computation of curvilinear coordinates
s=Calcdist1(X,Y,num,numtot);
ends=zeros(1,4);ends(1)=min(X);ends(2)=max(X);ends(3)=min(Y);ends(4)=max(Y);
[X,Y,~,num,numtot,~]=SplinePoints(deltas,lim,X,Y,s,num,numtot,lung);
s=Calcdist1(X,Y,num,numtot);


% Computation of initial curvature
[C,lim,ind]=Calccurv2(deltas,X,Y,s,num,numtot,primocut);

tempo=0;

% Initialization savings
if saveflag == 1
    params=['tau0=' num2str(tau0) ' ds=' num2str(ds) ' beta=' num2str(beta) ' deltat=' num2str(deltat) ' D0=' num2str(D0) ' T0=' num2str(T0) 'E' num2str(E)];
    name1=[namefile 'a'];name2=[namefile 'b'];name3=[namefile 'c'];
    name4=[namefile 'd'];name5=[namefile 'e'];name6=[namefile 'f'];
    nome7bin=[namefile 'values_occ'];
    save(name1,'params'); save(name2,'params'); save(name3,'params');
    save(name4,'params'); save(name5,'params'); save(name6,'params');

    lengthmaxcutoff=5000; %Maximum cutoff length
    cutoff.X=NaN(1,lengthmaxcutoff+1);cutoff.Y=NaN(1,lengthmaxcutoff+1);
    cutoff.s=NaN(1,lengthmaxcutoff+1);cutoff.C=NaN(1,lengthmaxcutoff+1);
    save(name1,'-append','cutoff');save(name2,'-append','cutoff');
    save(name3,'-append','cutoff');save(name4,'-append','cutoff');
    save(name5,'-append','cutoff');save(name6,'-append','cutoff');
end
if saveflag == 1 || saveflag == 2
    if qflag==1
        name8=[namefile 'S-V_mean-q'];save(name8,'params');
        data.t=0;data.S=1;data.V_mean=0;data.lc=0;
        save(name8,'-append','data');
    elseif qflag==2
        name8=[namefile 'S-V_mean-q'];save(name8,'params');
        data.t=0;data.S=1;data.V_mean=0;data.lc=0;data.q=0;
        save(name8,'-append','data');
    end
end

if qflag == 1
    q = Q0;
end

% Time step control
ti=1; ti2=1; S_iniz=0;
Ls=s(num);Lx=X(num);


for time=time0:deltat:t
    tempo=tempo+1;
    if qflag == 2
        q = Qstoch(tempo);
        if limitflag == 1 %upper limit
            if q > 2*Q0
                q = 2*Q0;
            end
        end
    end

    [H,~]=Compute_H(q,d50,I0,b,Hmax,dH,Delta,bedload);
    ds=d50/H; theta=I0/(Delta*ds);
    [trasp,~,~,~,~,~,~,~]=resistance(1,bedload,theta,ds,Delta,d50);

    Xv=X;Yv=Y;Cv=C;sv=s;numv=num;
    %%%%%%%%% Compute of the excess bank velocity and shift points 
    [X,Y,C,s,num,numtot,ind,lim,S,V_mean,v] = Itcurv1(X,Y,C,s,modelflag,b,num,bedload,q,I0,Hmax,dH,d50,deltat,deltas,E,H0,Cf0,numtot,primocut,Delta,2,Lx,numeff,saveflag,time,name8,ncutoff,lapse,savings,mainFolder,qflag);
    if max(abs(v*deltat/b))>0.1
        deltathalf=deltat/2;
        [X,Y,C,s,num,numtot,ind,lim,S,V_mean,v] = Itcurv1(Xv,Yv,Cv,sv,modelflag,b,numv,bedload,q,I0,Hmax,dH,d50,deltathalf,deltas,E,H0,Cf0,numtot,primocut,Delta,2,Lx,numeff,saveflag,time,name8,ncutoff,lapse,savings,mainFolder,qflag);
        [X,Y,C,s,num,numtot,ind,lim,S,V_mean,v] = Itcurv1(X,Y,C,s,modelflag,b,num,bedload,q,I0,Hmax,dH,d50,deltat,deltas,E,H0,Cf0,numtot,primocut,Delta,2,Lx,numeff,saveflag,time,name8,ncutoff,lapse,savings,mainFolder,qflag);
    end
    ti=ti+1;

    % Cut-off check
    ends(1)=min(X);ends(2)=max(X);ends(3)=min(Y);ends(4)=max(Y);
    numrig=(floor((ends(4)-ends(3))/dim))+1;
    numcol=(floor((ends(2)-ends(1))/dim))+1;
    numcell=numrig*numcol;

    ptimax=round(sqrt(2)*dim/(deltas))+2;

    Xmin=100;nn=find(X>Xmin); lc=s(num)-s(nn(1));
    [X,Y,C,s,lim,IND0,num,primocut,ncutoff]=Control_Cutoff(C,b,ends,numcell,numcol,numrig,ptimax,deltas,dim,X,Y,s,num,primocut,time,lim,ncutoff,lc,saveflag,numtot,name1,name2,name3,name4,name5,name6);


    while (IND0==1)
        [C,lim,ind]=Calccurv2(deltas,X,Y,s,num,numtot,primocut);
        [X,Y,C,s,lim,IND0,num,primocut,ncutoff]=Control_Cutoff(C,b,ends,numcell,numcol,numrig,ptimax,deltas,dim,X,Y,s,num,primocut,time,lim,ncutoff,lc,saveflag,numtot,name1,name2,name3,name4,name5,name6);
    end

    if cutoffmaxflag == 0 %initial configuration
        if S > 1.1
            nome_conf=['conf_end_' namefile]; save(nome_conf,'X','Y');
            break
        end
    elseif cutoffmaxflag == 1 %short-term simulations
        if (ncutoff >= Ncutoff_max) && (S > 1.5)
            break
        end
    end

end


text=['N. cutoff=' num2str(ncutoff)];disp(text);
text2=['Time=' num2str(floor(time/365)) ' years'];disp(text2);
text3=['Game over ' namefile];disp(text3);


