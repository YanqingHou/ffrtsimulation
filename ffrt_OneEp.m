function [resffrt]=ffrt_OneEp(ep,Qa,Nsamp)%Qa,Qb,Qab,Ps
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

counts=struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],...
    'themus',[],'thePfix',[],'thePfs',[],'thePscons',[]);
% 1. FFRT simulation
% 1.1 Initialization for FFRT
mus=0.001:0.001:1;
mus=mus';
muslen=length(mus);

Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';
Pfslen=length(Pfs);

ffrtclnep=1;ffrtclnna=2;ffrtclnns=3;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;
ffrtclnmu=7;ffrtclnPfix=8;ffrtclnPscon=9; ffrtclnPftrue=10;

resffrt=zeros(na*Pfslen,10);


% 1.2 initialization of the struct array
for ns= 1:na
    k=na-ns+1;
    counts(ns).Psb=prod(2 * normcdf(1./(2*sqrt(D(k:end)))) -1 );
    counts(ns).sucnumILS=0;
   
    counts(ns).fixnumrt=zeros(muslen,1);
    counts(ns).sucnumrt=counts(ns).fixnumrt;
    counts(ns).failnumrt=counts(ns).fixnumrt;

end


% 1.3 Simulate all the samples and do LAMBDA, ratio test.
for i=1:Nsamp
        zhat=mvnrnd(ztrue,Qzhat,1)';
    for ns=1:na
        k=na-ns+1; 
        [zpar,sqnorm] = ssearch(zhat(k:end),L(k:end,k:end),D(k:end),ncands);
        zpar=zpar(:,1);
        
        sucflag=sum(~(ztrue(k:end)==zpar))==0;
        counts(ns).sucnumILS = counts(ns).sucnumILS+sucflag;         % Correctly fixed

        counts(ns).ratio=sqnorm(1)/sqnorm(2);
        
        rtpass=counts(ns).ratio<mus;                % 100 different mu
        counts(ns).fixnumrt=counts(ns).fixnumrt+rtpass;             % accepted by FFRT
        counts(ns).sucnumrt=counts(ns).sucnumrt + double(rtpass&sucflag);   % success rate after FFRT: (Accepted by FFRT & correct) / N
        counts(ns).failnumrt=counts(ns).fixnumrt-counts(ns).sucnumrt;% failure rate after FFRT (Fix-Fix & correct)
    end
 
end

% parameter for FFRT
for ns=1:na
    
    counts(ns).sucnumrt=counts(ns).sucnumrt./counts(ns).fixnumrt; % Pscon
    counts(ns).failnumrt=counts(ns).failnumrt/Nsamp; % Pf
    counts(ns).fixnumrt=counts(ns).fixnumrt/Nsamp; % Pfix
    % if the fix rate is 0, set failure rate to 1 and conditaional success
    % rate to 0.
    onesidx=counts(ns).fixnumrt==0; 
    counts(ns).failnumrt(onesidx)=1;
    counts(ns).sucnumrt(onesidx)=0;
    for pfi=1:Pfslen
        
        [rowtmp, ~, mupf]=find(counts(ns).failnumrt<=Pfs(pfi),1,'last'); % find the critical value with a given failure rate
        if isempty(mupf)
%            mupf=0; 

           counts(ns).themus(pfi)=0;%mupf; mu=0 means reject any candidate
           counts(ns).thePfix(pfi)=0;%counts(ns).fixnumrt(row);
           counts(ns).thePscons(pfi)=0;%counts(ns).sucnumrt(row);
           counts(ns).thePfs(pfi)=0;%counts(ns).failnumrt(row);

        else
            counts(ns).themus(pfi)=mus(rowtmp);
            counts(ns).thePfix(pfi)=counts(ns).fixnumrt(rowtmp);
            counts(ns).thePscons(pfi)=counts(ns).sucnumrt(rowtmp);
            counts(ns).thePfs(pfi)=counts(ns).failnumrt(rowtmp);
                    
        end
         
         
    end
% counts=struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],...
%     'themus',[],'thePfix',[],'thePfs',[],'thePscons',[]);
%     resffrt(ns,:)=[ep,counts(ns).Pftruert(pfi),];

        indffrt=(ns-1)*Pfslen+1:ns*Pfslen;
        resffrt(indffrt,ffrtclnep)=ep; 
        resffrt(indffrt,ffrtclnna)=na; 

        resffrt(indffrt,ffrtclnns)=ns;
        resffrt(indffrt,ffrtclnPsb)=counts(ns).Psb;
        resffrt(indffrt,ffrtclnPsILS)=counts(ns).sucnumILS/Nsamp;
        resffrt(indffrt,ffrtclnPf_req)=Pfs;
        resffrt(indffrt,ffrtclnmu)=counts(ns).themus;
        resffrt(indffrt,ffrtclnPfix)=counts(ns).thePfix;
        resffrt(indffrt,ffrtclnPscon)=counts(ns).thePscons;
        
        resffrt(indffrt,ffrtclnPftrue)=counts(ns).thePfs;
      
end
