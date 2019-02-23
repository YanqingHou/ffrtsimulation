function [resffrt]=valid_fit_ffrt_OneEp(ep,Qa,Nsamp)%Qa,Qb,Qab,Ps
% function [reseprt01, reseprt001, resepdt01, resepdt001]=muFixRateOneEp(ep,Qa,Qb,Qab,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
%
% TSRConeep(ep,Qa,Qb,Qab,Psb,Nsamp)%Qa,Qb,Qab,Ps
% input: Qa, Qb, Qab, Ps

% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
na      = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
ncands=2;
[Qzhat,~,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');

% 1. FFRT simulation
% 1.1 Initialization for FFRT
% mus=0.001:0.001:1;
% mus=mus';
% muslen=length(mus);
Pfs=[0.0005, 0.0006, 0.0007, 0.0008, 0.0009,...
    0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01];
Pfslen=length(Pfs);

ffrtclnep=1;ffrtclnna=2;ffrtclnns=3;ffrtclnPsb=4;ffrtclnPf_req=5;
ffrtclnmtd=6;ffrtclnmu=7;ffrtclnPfix=8;ffrtclnPscon=9; ffrtclnPftrue=10; ffrtclnfalsealarm=11;

mtdnum=6;
% methods: 1--mu=1; 2--mu=1.5; 3--mu=2; 4--mu=3; 5--mu=sandra table; 6--mu=fitting funcs;
resffrt=zeros(na*Pfslen*mtdnum,11);

counts=struct('Psb',[],'ratio',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsenumrt',[],...
    'mus',[],'Pfix',[],'Pfs',[],'Pscons',[],'falsealarm',[]);

% 1.2 initialization of the struct array
for ns= 1:na
    k=na-ns+1;
    counts(ns).Psb=prod(2 * normcdf(1./(2*sqrt(D(k:end)))) -1 );
    %     counts(ns).sucnumILS=0;
    
    counts(ns).fixnumrt=zeros(mtdnum,Pfslen);
    counts(ns).sucnumrt=counts(ns).fixnumrt;
    counts(ns).failnumrt=counts(ns).fixnumrt;
    counts(ns).falsenumrt=counts(ns).fixnumrt;
    counts(ns).mus=counts(ns).fixnumrt;
    for ipf=1:Pfslen
        Pfsv=0.001; if Pfs(ipf)>0.001, Pfsv=0.01; end
        svmu=svcpmu(1-counts(ns).Psb,Pfsv,ns);%cpsandramu();
        Fitmu=cpmu1(1-counts(ns).Psb,Pfs(ipf),ns);%7;%cpmu1();
        counts(ns).mus(:,ipf)=[1,1/1.5,1/2,1/3,svmu,Fitmu]';
    end
end


% 1.3 Simulate all the samples and do LAMBDA, ratio test.
for i=1:Nsamp
    zhat=mvnrnd(ztrue,Qzhat,1)';
    for ns=1:na
        k=na-ns+1;
        [zpar,sqnorm] = ssearch(zhat(k:end),L(k:end,k:end),D(k:end),ncands);
        zpar=zpar(:,1);
        
        sucflag=sum(~(ztrue(k:end)==zpar))==0;
        %         counts(ns).sucnumILS = counts(ns).sucnumILS+sucflag;         % Correctly fixed
        
        counts(ns).ratio=sqnorm(1)/sqnorm(2);
        
        rtpass=counts(ns).ratio<=counts(ns).mus;                % 6 different mu
        counts(ns).fixnumrt=counts(ns).fixnumrt+rtpass;             % accepted by FFRT
        counts(ns).sucnumrt=counts(ns).sucnumrt + double(rtpass&sucflag);   % success rate after FFRT: (Accepted by FFRT & correct) / N
        counts(ns).failnumrt=counts(ns).fixnumrt-counts(ns).sucnumrt;% failure rate after FFRT (Fix-Fix & correct)
        counts(ns).falsenumrt=counts(ns).falsenumrt+double((~rtpass)&sucflag);% false alarm after FFRT (Fix-Fix & correct)
    end
    
end

% parameter for FFRT
for ns=1:na
    
    counts(ns).sucnumrt=counts(ns).sucnumrt./counts(ns).fixnumrt; % Pscon
    counts(ns).failnumrt=counts(ns).failnumrt/Nsamp; % Pf
    counts(ns).fixnumrt=counts(ns).fixnumrt/Nsamp; % Pfix
    counts(ns).falsenumrt=counts(ns).falsenumrt/Nsamp;
    idxnan=isnan(counts(ns).sucnumrt);
    counts(ns).sucnumrt(idxnan)=0;
    
    for ipf=1:Pfslen
        indffrt=(ns-1)*(Pfslen*mtdnum)+(ipf-1)*mtdnum+1:(ns-1)*(Pfslen*mtdnum)+ipf*mtdnum;
        
        
        resffrt(indffrt,ffrtclnep)=ep;
        resffrt(indffrt,ffrtclnna)=na;
        
        resffrt(indffrt,ffrtclnns)=ns;
        resffrt(indffrt,ffrtclnPsb)=counts(ns).Psb;
        %         resffrt(indffrt,ffrtclnPsILS)=counts(ns).sucnumILS/Nsamp;
        resffrt(indffrt,ffrtclnPf_req)=Pfs(ipf);
        resffrt(indffrt,ffrtclnmtd)=[1,2,3,4,5,6]';
        
        resffrt(indffrt,ffrtclnmu)=counts(ns).mus(:,ipf);
        resffrt(indffrt,ffrtclnPfix)=counts(ns).fixnumrt(:,ipf);
        resffrt(indffrt,ffrtclnPscon)=counts(ns).sucnumrt(:,ipf);
        
        resffrt(indffrt,ffrtclnPftrue)=counts(ns).failnumrt(:,ipf);
        resffrt(indffrt,ffrtclnfalsealarm)=counts(ns).falsenumrt(:,ipf);
    end
    
end
