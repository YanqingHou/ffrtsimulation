function [reseprt01, reseprt001, resepdt01, resepdt001]=FixRateOneEp(ep,Qa,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
% clnep=1;clnns=2;clnPs=3;clnPfix=4;clnPsb=5;clnPf=6;
ncands=2;
na  = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
[Qzhat,~,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');
if P0==0
    n0=na;
else
    n0=Ps2ns(Psb,P0,na,D);
end

if n0==0 
    reseprt01=[ep,zeros(1,5)];
    reseprt001=reseprt01;
    resepdt01=reseprt01;
    resepdt001=reseprt01;
    return ; 
end
Pss=zeros(n0,1);

fixnum01=Pss;
sucnum01=Pss;
failnum01=Pss;

fixnum001=Pss;
sucnum001=Pss;
failnum001=Pss;

fixnumdt01=Pss;
sucnumdt01=Pss;
failnumdt01=Pss;


fixnumdt001=Pss;
sucnumdt001=Pss;
failnumdt001=Pss;

mu01s=Pss;
mu001s=Pss;
% sucflag=0;
for ns=1:n0
    k1=na-ns+1;
    Pss(ns)= prod (2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
    if Pss(ns)>=1-0.01
        mu01s(ns)=1;
    else
        mu01s(ns)=ratioinv(0.01,1-Pss(ns),ns); 
    end
    
    if Pss(ns)>=1-0.001,
        mu001s(ns)=1;
    else
        mu001s(ns)=ratioinv(0.001,1-Pss(ns),ns); 
    end
end

for i=1:Nsamp
        zhat=mvnrnd(ztrue,Qzhat,1)';
    for ns=1:n0
        k1=na-ns+1; %zp=zhat(na-ns+1:na);
        [zpar1,sqnorm] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);
        zpar1=zpar1(:,1);
        sucflag=sum(~(ztrue(k1:end)==zpar1))==0;
%  FFRT with Pf=0.01
       if sqnorm(1)/sqnorm(2)<mu01s(ns)%s FFRT in the first step accept
               fixnum01(ns)=fixnum01(ns)+1;
            if sucflag
               sucnum01(ns)=sucnum01(ns)+1;
            else
               failnum01(ns)=failnum01(ns)+1;    
            end
       end
        
%         FFRT with Pf=0.001
        if sqnorm(1)/sqnorm(2)<mu001s(ns)%s FFRT in the first step accept
               fixnum001(ns)=fixnum001(ns)+1;
            if sucflag
               sucnum001(ns)=sucnum001(ns)+1;
            else
               failnum001(ns)=failnum001(ns)+1;    
            end
        end 
        % FFDT with Pf=0.01
        if sqnorm(2)-sqnorm(1)>=ffdtinv(0.01,Pss(ns))
               fixnumdt01(ns)=fixnumdt01(ns)+1;
            if sucflag
               sucnumdt01(ns)=sucnumdt01(ns)+1;
            else
               failnumdt01(ns)=failnumdt01(ns)+1;    
            end 
        end
        
       % FFDT with Pf=0.001
         if  sqnorm(2)-sqnorm(1)>=ffdtinv(0.001,Pss(ns))
               fixnumdt001(ns)=fixnumdt001(ns)+1;
            if sucflag
               sucnumdt001(ns)=sucnumdt001(ns)+1;
            else
               failnumdt001(ns)=failnumdt001(ns)+1;    
            end 
        end
    end
end
nss=1:n0;nss=nss';
Pfix=fixnum01/Nsamp; Pscon=sucnum01./fixnum01; Pf=failnum01/Nsamp; nanidx=isnan(Pscon); Pscon(nanidx)=0; reseprt01=[ep*ones(n0,1),nss,Pss,Pfix,Pscon,Pf];
Pfix=fixnum001/Nsamp; Pscon=sucnum001./fixnum001; Pf=failnum001/Nsamp; nanidx=isnan(Pscon); Pscon(nanidx)=0; reseprt001=[ep*ones(n0,1),nss,Pss,Pfix,Pscon,Pf];
Pfix=fixnumdt01/Nsamp; Pscon=sucnumdt01./fixnumdt01; Pf=failnumdt01/Nsamp; nanidx=isnan(Pscon); Pscon(nanidx)=0; resepdt01=[ep*ones(n0,1),nss,Pss,Pfix,Pscon,Pf];
Pfix=fixnumdt001/Nsamp; Pscon=sucnumdt001./fixnumdt001; Pf=failnumdt001/Nsamp; nanidx=isnan(Pscon); Pscon(nanidx)=0; resepdt001=[ep*ones(n0,1),nss,Pss,Pfix,Pscon,Pf];