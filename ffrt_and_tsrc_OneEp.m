function [resffrt,restsrc]=ffrt_and_tsrc_OneEp(ep,Qa,Qb,Qab,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
% function [reseprt01, reseprt001, resepdt01, resepdt001]=muFixRateOneEp(ep,Qa,Qb,Qab,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
% 
% TSRConeep(ep,Qa,Qb,Qab,Psb,Nsamp)%Qa,Qb,Qab,Ps
% input: Qa, Qb, Qab, Ps

% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
na      = size(Qa,1);
nb      = size(Qb,1);
Qx=[Qa,Qab;
    Qab', Qb];
Qx      = (tril(Qx,0)+tril(Qx,-1)');

atrue=randi([-200 200],na,1);% generate random integers
btrue=zeros(nb,1);
xtrue=[atrue;btrue];
% clnep=1;clnns=2;clnPs=3;clnPfix=4;clnPsb=5;clnPf=6;
ncands=2;
% na  = size(Qa,1);
% atrue=randi([-200 200],na,1);% generate random integers
[Qzhat,Z,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');


% na=na;
% if P0>0
%     n0=Ps2ns(Psb,P0,na,D);
% end
counts=struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],...
    'themus',[],'thePfix',[],'thePfs',[],'thePscons',[]);
%            counts(ns).muPfs(pfi)=0;%mupf; mu=0 means reject any candidate
%            counts(ns).Pfixrts(pfi)=0;%counts(ns).fixnumrt(row);
%            counts(ns).Psconrts(pfi)=0;%counts(ns).sucnumrt(row);
%            counts(ns).muavail(pfi,:)=counts(ns).fltavail; %counts(ns).availsrt(muslen,:);% should be the availability of the float solution
%            counts(ns).Pftruert(pfi)=0;%counts(ns).failnumrt(row);
% 1. FFRT simulation
% 1.1 Initialization for FFRT
mus=0.01:0.01:1;
mus=mus';
muslen=length(mus);

Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';
Pfslen=length(Pfs);

resffrt=zeros(na*Pfslen,11);


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

%all the possible result simulated. In the next, we need to find the best
%parameter for FFRT and TSRC

ffrtclnep=1;ffrtclnna=2;ffrtclnns=3;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;
ffrtclngain=7;ffrtclnmu=8;ffrtclnPfix=9;ffrtclnPscon=10;
ffrtclnPftrue=11;

% 
% Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';

rows=zeros(na,Pfslen);
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
    %row=zeros(Pfslen,1);
    for pfi=1:Pfslen
        
        [rowtmp, ~, mupf]=find(counts(ns).failnumrt<=Pfs(pfi),1,'last'); % find the critical value with a given failure rate
        if isempty(mupf)
%            mupf=0; 

           counts(ns).themus(pfi)=0;%mupf; mu=0 means reject any candidate
           counts(ns).thePfix(pfi)=0;%counts(ns).fixnumrt(row);
           counts(ns).thePscons(pfi)=0;%counts(ns).sucnumrt(row);
           counts(ns).thePfs(pfi)=0;%counts(ns).failnumrt(row);

        else
%             rows(ns,pfi)=rowtmp;
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

%Pfpairs=[0.001, 0.0005;
 %        0.001, 0.0006;
  %       0.001, 0.0007;
   %      0.001, 0.0008;
    %     0.001, 0.0009;
     %    0.01, 0.001;
      %   0.01, 0.002;
       %  0.01, 0.003;
        % 0.01, 0.004;
        % 0.01, 0.005;
        % 0.01, 0.006;
        % 0.01, 0.007;
        % 0.01, 0.008;
        % 0.01, 0.009];
        Pfpairslen=size(Pfpairs,1);
        indrestsrc=0;
    for iPfpair=1:Pfpairslen
        Pfttl=Pfpairs(iPfpair,1);
        if Pfttl==0.001, indpf1=6; else indpf1=15; end
        Pf1=Pfpairs(iPfpair,2);
        ns2=Ps2ns(Psb,1-Pf1,na,D);
        for ns=ns2+1:na
%             indtsrc=Pfpairslen*()
            indrestsrc=indrestsrc+1;
            Pfix1=counts(ns).Pfixrts(iPfpair);
            if Pfix1==0
% 
%                counts(ns).muPfs(pfi)=mus(rowtmp);
%             counts(ns).Pfixrts(pfi)=counts(ns).fixnumrt(rowtmp);
%             counts(ns).Psconrts(pfi)=counts(ns).sucnumrt(rowtmp);
%             counts(ns).muavail(pfi,:)=counts(ns).availsrt(rowtmp,:);
%             counts(ns).Pftruert(pfi)=counts(ns).failnumrt(rowtmp);
% %                       mu2=1;
%                         Pf2true=counts(ns).mufailnum2(iPfpair,muslen);
%                         muavails2=avails2{rows(ns,iPfpair),iPfpair,ns}(muslen,:);
%                         Pfix2=0;
%                         Pscon2=0;
                        
%               Pf2=Pfttl;
              if ns2==0
                  mu2=1;
                  Pf2true=0;
                  Pfix2=0;
                  Pscon2=0;
                  muavails2=comparevect(beflt,errs);
              else
                  mu2=counts(ns2).muPfs(indpf1);
                  Pf2true=counts(ns2).Pftruert(iPfpair);
                  Pfix2=counts(ns2).Pfixrts(iPfpair);
                  Pscon2=counts(ns2).Psconrts(iPfpair);
                  muavails2=counts(ns2).muavail(iPfpair,:);
              end
          elseif Pfix1<1
              if ns2==0
                  mu2=1;
                  Pf2true=0;
                  Pfix2=0;
                  Pscon2=0;
                  muavails2=comparevect(beflt,errs);                  
              else
                    Pf2=(1-Pf1)/(1-Pfix1);   
                    if Pf2<1

                        counts(ns).musucnum2(iPfpair,:)=counts(ns).sucnum2(:,rows(ns,iPfpair),iPfpair)./counts(ns).fixnum2(:,rows(ns,iPfpair),iPfpair); nanidx=isnan(counts(ns).musucnum2(iPfpair,:)); counts(ns).musucnum2(iPfpair,nanidx)=0; 
                        counts(ns).mufixnum2(iPfpair,:)=counts(ns).fixnum2(:,rows(ns,iPfpair),iPfpair)./(Nsamp*(1-counts(ns).Pfixrts(iPfpair))); nanidx=isnan(counts(ns).mufixnum2(iPfpair,:)); counts(ns).mufixnum2(iPfpair,nanidx)=0; 
                        counts(ns).mufailnum2(iPfpair,:)=counts(ns).failnum2(:,rows(ns,iPfpair),iPfpair)./counts(ns).fixnum2(:,rows(ns,iPfpair),iPfpair); nanidx=isnan(counts(ns).mufailnum2(iPfpair,:)); counts(ns).mufailnum2(iPfpair,nanidx)=1;

                        %Pfttlemp=counts(ns).failnum2(:,rows(ns,iPfpair),iPfpair)/Nsamp+Pf1;
                       [rowmu2,~, Pf2true]=find(counts(ns).mufailnum2(iPfpair,:)<=Pf2,1,'last');
                        mu2=mus(rowmu2);
                        muavails2=avails2{rows(ns,iPfpair),iPfpair,ns}(rowmu2,:);
                        Pfix2=counts(ns).mufixnum2(iPfpair,rowmu2);
                        Pscon2=counts(ns).musucnum2(iPfpair,rowmu2);
                    else
                        mu2=1;
                        Pf2true=counts(ns).mufailnum2(iPfpair,muslen);
                        muavails2=zeros(1,errslen);%avails2{rows(ns,iPfpair),iPfpair,ns}(muslen,:);
                        Pfix2=0;
                        Pscon2=0;
                    end
              end
            else 
                mu2=1;
                Pf2true=0;
                muavails2=zeros(1,errslen);%avails2{rows(ns,iPfpair),iPfpair,ns}(muslen,:);
                Pfix2=0;
                Pscon2=0;

            end
            
         availstsrc=counts(ns).muavail(iPfpair,:) + muavails2./Nsamp;
            
         mu1=counts(ns).muPfs(indpf1);
         Pfix1=counts(ns).Pfixrts(indpf1);
         Pscon1=counts(ns).Psconrts(indpf1);
         Pftrue1=counts(ns).Pftruert(indpf1);
         if ns2==0
            Pssub2=0;
            g2=1;
         else
            Pssub2=counts(ns2).Psb;
            g2= counts(ns2).gain;  
         end
         restsrc(indrestsrc,:)=[ep, na, Psb, ns, counts(ns).Psb, counts(ns).gain, Pf1, mu1, Pfix1, Pscon1, Pftrue1, ns2, Pssub2, g2, mu2, Pfix2, Pscon2, Pf2true, availstsrc ];   

%          restsrc=[restsrc; ep, na, Psb, ns, counts(ns).Psb, counts(ns).gain, Pf1, mu1, Pfix1, Pscon1, Pftrue1, ns2, counts(ns2).Psb, counts(ns2).gain, mu2, Pfix2, Pscon2, Pf2true, availstsrc ];   
         %[rowmu21,~,Pf2true]=stdbhat,

        end
    end
