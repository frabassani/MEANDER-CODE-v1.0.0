function SaveCutoff(Xt,Yt,st,Ct,name1,name2,name3,name4,name5,name6,time)

switch (floor(time/25))
      case 0
            pp=load(name1);
            firstrow=pp.cutoff.X;
            lengthmaxcutoff=length(firstrow)-1;
            cutoffX=[pp.cutoff.X;[Xt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffY=[pp.cutoff.Y;[Yt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffs=[pp.cutoff.s;[st NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffC=[pp.cutoff.C;[Ct NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoff.X=cutoffX;cutoff.Y=cutoffY;cutoff.s=cutoffs;cutoff.C=cutoffC;
            save(name1,'-append','cutoff');
      case 1
            pp=load(name2);
            firstrow=pp.cutoff.X;
            lengthmaxcutoff=length(firstrow)-1;
            cutoffX=[pp.cutoff.X;[Xt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffY=[pp.cutoff.Y;[Yt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffs=[pp.cutoff.s;[st NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffC=[pp.cutoff.C;[Ct NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoff.X=cutoffX;cutoff.Y=cutoffY;cutoff.s=cutoffs;cutoff.C=cutoffC;
            save(name2,'-append','cutoff');    
      case 2
            pp=load(name3);
            firstrow=pp.cutoff.X;
            lengthmaxcutoff=length(firstrow)-1;
            cutoffX=[pp.cutoff.X;[Xt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffY=[pp.cutoff.Y;[Yt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffs=[pp.cutoff.s;[st NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffC=[pp.cutoff.C;[Ct NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoff.X=cutoffX;cutoff.Y=cutoffY;cutoff.s=cutoffs;cutoff.C=cutoffC;
            save(name3,'-append','cutoff');    
      case 3
            pp=load(name4);
            firstrow=pp.cutoff.X;
            lengthmaxcutoff=length(firstrow)-1;
            cutoffX=[pp.cutoff.X;[Xt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffY=[pp.cutoff.Y;[Yt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffs=[pp.cutoff.s;[st NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffC=[pp.cutoff.C;[Ct NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoff.X=cutoffX;cutoff.Y=cutoffY;cutoff.s=cutoffs;cutoff.C=cutoffC;
            save(name4,'-append','cutoff');
      case 4
            pp=load(name5);
            firstrow=pp.cutoff.X;
            lengthmaxcutoff=length(firstrow)-1;
            cutoffX=[pp.cutoff.X;[Xt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffY=[pp.cutoff.Y;[Yt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffs=[pp.cutoff.s;[st NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffC=[pp.cutoff.C;[Ct NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoff.X=cutoffX;cutoff.Y=cutoffY;cutoff.s=cutoffs;cutoff.C=cutoffC;
            save(name5,'-append','cutoff');
      otherwise
            pp=load(name6);
            firstrow=pp.cutoff.X;
            lengthmaxcutoff=length(firstrow)-1;
            cutoffX=[pp.cutoff.X;[Xt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffY=[pp.cutoff.Y;[Yt NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffs=[pp.cutoff.s;[st NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoffC=[pp.cutoff.C;[Ct NaN(1,lengthmaxcutoff-length(Xt)+1)]];
            cutoff.X=cutoffX;cutoff.Y=cutoffY;cutoff.s=cutoffs;cutoff.C=cutoffC;
            save(name6,'-append','cutoff');
end 

end %function