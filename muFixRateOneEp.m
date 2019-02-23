function [resffrt,restsrc]=muFixRateOneEp(ep,Qa,Qb,Qab,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
% function [reseprt01, reseprt001, resepdt01, resepdt001]=muFixRateOneEp(ep,Qa,Qb,Qab,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
% 
% TSRConeep(ep,Qa,Qb,Qab,Psb,Nsamp)%Qa,Qb,Qab,Ps
% input: Qa, Qb, Qab, Ps
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

n0=na;
if P0>0
    n0=Ps2ns(Psb,P0,na,D);
end


counts=struct('ns',[],'Psb',[],'Kcoef1',[],'gain',[],'errs',[],'availsrt',[],...
    'ratio',[],'sucnumILS',[],'fltavail',[],...
    'fixnumrt',[],'sucnumrt',[],'failnumrt',[],...
    'Pfs',[],'muPfs',[],...
    'Pfixrts',[],'Psconrts',[],...
    'muavail',[],'Pftruert',[],...
    'fixnum2',[],'sucnum2',[],'failnum2',[],...
    'mufixnum2',[],'musucnum2',[],'mufailnum2',[]);
% Pss=zeros(n0,1);
mus=0.01:0.01:1;
mus=mus';
muslen=length(mus);
errs=[0.01:0.01:0.1,0.2:0.1:0.5]';
errslen=length(errs);
Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';
Pfslen=length(Pfs);
% avails=zeros(muslen,errslen);

if n0==0 
    resffrt=zeros(1,11);
    resffrt(1)=ep;

    restsrc=zeros(1,18+errslen);
    restsrc(1)=ep;    
    return; 
end
resffrt=zeros(n0*Pfslen,11);
% restsrc=zeros(n0*Pfslen*muslen,19+errslen+1);% one more row for float solution
restsrc=zeros(n0*Pfslen,18+errslen);% one more row for float solution
restsrc(:,1)=ep;

avails2=cell(muslen,Pfslen,n0);
Ps0=zeros(1,Pfslen);
ns2=Ps0;
for i=1:Pfslen
Ps0(i)=1-Pfs(i);
ns2(i)=Ps2ns(Psb,Ps0(i),na,D);
end

% initialization of the struct array
for ns= 1:n0
    k=na-ns+1;
    counts(ns).ns=ns;
    counts(ns).Psb=prod(2 * normcdf(1./(2*sqrt(D(k:end)))) -1 );
    counts(ns).sucnumILS=0;
    counts(ns).fltavail = zeros(1,errslen); 
    
    Qzpar1 = Qzhat(k:end,k:end); Zpar1 = Z(:,k:end); Qbzpar1=Qab'*Zpar1;
    Kcoef1=Qbzpar1/Qzpar1;  
    Qbp1=Qb-Kcoef1*Qbzpar1';
    Rs1=sqrt(det(Qb(1:3,1:3))/det(Qbp1(1:3,1:3)))^(1/3);%large than 1
    counts(ns).gain=Rs1;
    counts(ns).Kcoef1=Kcoef1;
    
    counts(ns).fixnumrt=zeros(muslen,1);
    counts(ns).sucnumrt=counts(ns).fixnumrt;
    counts(ns).failnumrt=counts(ns).fixnumrt;
  
    %counts(ns).errs=errs;  % for saving in the output file
    counts(ns).availsrt=zeros(muslen,errslen); % for each mu and each required err, there is an availability
    counts(ns).muavail=zeros(Pfslen,errslen);  % There are 11 Pfs, for each required Pf, find the mu value; and for this mu value, calculate the availability of 20 errors. This matrix is a part of the matrix availsrt.


    %%%% The result for TSRC, needs to be confirmed. 
    counts(ns).Pfs=Pfs;
    counts(ns).fixnum2=zeros(muslen,muslen,Pfslen); % for each subset, apply different Pfs, and for each Pf, ....
    counts(ns).sucnum2=zeros(muslen,muslen,Pfslen);
    counts(ns).failnum2=zeros(muslen,muslen,Pfslen);
    
    counts(ns).musucnum2=zeros(Pfslen,muslen);
    counts(ns).mufixnum2=zeros(Pfslen,muslen);
    counts(ns).mufailnum2=zeros(Pfslen,muslen);

    counts(ns).muPfs=zeros(Pfslen,1);
    counts(ns).Pfixrts=zeros(Pfslen,1);
    counts(ns).Psconrts=zeros(Pfslen,1);
    counts(ns).Pftruert=zeros(Pfslen,1);


    for i=1:Pfslen
        for j=1:muslen
            avails2{j,i,ns}=zeros(muslen,errslen);
        end
    end

end



for i=1:Nsamp
        xhat=mvnrnd(xtrue,Qx,1)';
        ahat=xhat(1:na); bhat=xhat(na+1:na+nb);
        zhat=Z'*ahat;
        bert=norm(bhat(1:3))*ones(muslen,n0);   
        beflt=norm(bhat(1:3));
%         berttsrc=bert;   

%         zhat=mvnrnd(ztrue,Qzhat,1)';
        sucflags=zeros(n0,1);
    for ns=1:n0
        k=na-ns+1; 
        [zpar,sqnorm] = ssearch(zhat(k:end),L(k:end,k:end),D(k:end),ncands);
        zpar=zpar(:,1);
        bchk=bhat-counts(ns).Kcoef1*(zhat(k:end)-zpar);
        
        sucflag=sum(~(ztrue(k:end)==zpar))==0;
        counts(ns).sucnumILS = counts(ns).sucnumILS+sucflag;         % Correctly fixed
        sucflags(ns)=sucflag;
        counts(ns).fltavail = counts(ns).fltavail+double(bert(1)<=errs');
        counts(ns).ratio=sqnorm(1)/sqnorm(2);
        
        rtpass=counts(ns).ratio<mus;                % 100 different mu
        counts(ns).fixnumrt=counts(ns).fixnumrt+rtpass;             % accepted by FFRT
        counts(ns).sucnumrt=counts(ns).sucnumrt + double(rtpass&sucflag);   % success rate after FFRT: (Accepted by FFRT & correct) / N
        counts(ns).failnumrt=counts(ns).fixnumrt-counts(ns).sucnumrt;% failure rate after FFRT (Fix-Fix & correct)
        bert(rtpass,ns)=norm(bchk(1:3));                                % if FFRT passed, baseline estimation is updated, and the rms is updated.
        resmat=comparevect(bert(:,ns),errs);                    % compare the rms with required errors for each mu in FFRT, and this result in a matrix.
        counts(ns).availsrt=counts(ns).availsrt+resmat;                 % In the end, it needs to be devided by N---- the total sample number.
    end
    % Do TSRC here
   for j=1:Pfslen % for every Pf... wait, is it neccessary?
       if ns2(j)==0
           ;
       else
           
       rtpass2=counts(ns2(j)).ratio<=mus; % the ratio test of the second step
        for ns=ns2(j)+1:n0
            rtpass=counts(ns).ratio<=mus;                % 100 different mu
            rtreject=~rtpass';

            counts(ns).fixnum2(:,:,j)=counts(ns).fixnum2(:,:,j)+kron(rtpass2,rtreject); % the first index is for mus in step2, the second index is for mus in first step, the third index is for different Pf                                       %
            counts(ns).sucnum2(:,:,j)=counts(ns).sucnum2(:,:,j) + double( kron(rtpass2,rtreject)&sucflags(ns2(j)));   % success rate after FFRT: (Accepted by FFRT & correct) / N
            counts(ns).failnum2(:,:,j)=counts(ns).fixnum2(:,:,j)-counts(ns).sucnum2(:,:,j); 
                    

            for k=1:muslen  % for each mu in the first step of TSRC
                if rtreject(k)==1
                   resmat=comparevect(bert(:,ns2(j)),errs);                    % compare the rms with required errors for each mu in FFRT, and this result in a matrix.
                   avails2{k,j,ns}=avails2{k,j,ns}+resmat;% the subset in the first step is ns, the subset in the second step is ns2(j), the mu in the first step ratio test is mus(k); each element in avails2 is a matrix with muslen*errslen dimension. float �� fixed solution Ҫ�ֿ���
                end
            end

        end
        
   end
   end
end

%all the possible result simulated. In the next, we need to find the best
%parameter for FFRT and TSRC

ffrtclnep=1;ffrtclnna=2;ffrtclnns=3;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;
ffrtclngain=7;ffrtclnmu=8;ffrtclnPfix=9;ffrtclnPscon=10;
ffrtclnPftrue=11;

% Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';

Pfpairs=[0.001, 0.0005;
         0.001, 0.0006;
         0.001, 0.0007;
         0.001, 0.0008;
         0.001, 0.0009;
         0.01, 0.001;
         0.01, 0.002;
         0.01, 0.003;
         0.01, 0.004;
         0.01, 0.005;
         0.01, 0.006;
         0.01, 0.007;
         0.01, 0.008;
         0.01, 0.009;
         0.02,  0.01];

% 
% Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';

rows=zeros(n0,Pfslen);
% parameter for FFRT
for ns=1:n0
    
    counts(ns).sucnumrt=counts(ns).sucnumrt./counts(ns).fixnumrt; % Pscon
    counts(ns).failnumrt=counts(ns).failnumrt/Nsamp; % Pf
    counts(ns).fixnumrt=counts(ns).fixnumrt/Nsamp; % Pfix
    counts(ns).availsrt=counts(ns).availsrt/Nsamp; % avails in the first step of TSRC
    nanidx=isnan(counts(ns).sucnumrt);  counts(ns).sucnumrt(nanidx)=0;
    onesidx=counts(ns).fixnumrt==0; counts(ns).failnumrt(onesidx)=1;
    %row=zeros(Pfslen,1);
    for pfi=1:Pfslen
        
        [rowtmp, ~, mupf]=find(counts(ns).failnumrt<=Pfs(pfi),1,'last'); % find the critical value with a given failure rate
        if isempty(mupf)
%            mupf=0; 

           counts(ns).muPfs(pfi)=0;%mupf; mu=0 means reject any candidate
           counts(ns).Pfixrts(pfi)=0;%counts(ns).fixnumrt(row);
           counts(ns).Psconrts(pfi)=0;%counts(ns).sucnumrt(row);
           counts(ns).muavail(pfi,:)=counts(ns).fltavail; %counts(ns).availsrt(muslen,:);% should be the availability of the float solution
           counts(ns).Pftruert(pfi)=0;%counts(ns).failnumrt(row);
       
           counts(ns).musucnum2(pfi,:)=zeros(1,muslen);%for TSRC
           counts(ns).mufixnum2(pfi,:)=zeros(1,muslen);%for TSRC
           counts(ns).mufailnum2(pfi,:)=zeros(1,muslen);%for TSRC
        else
            rows(ns,pfi)=rowtmp;
            counts(ns).muPfs(pfi)=mus(rowtmp);
            counts(ns).Pfixrts(pfi)=counts(ns).fixnumrt(rowtmp);
            counts(ns).Psconrts(pfi)=counts(ns).sucnumrt(rowtmp);
            counts(ns).muavail(pfi,:)=counts(ns).availsrt(rowtmp,:);
            counts(ns).Pftruert(pfi)=counts(ns).failnumrt(rowtmp);
                    
        end
         
         
    end

%     resffrt(ns,:)=[ep,counts(ns).Pftruert(pfi),];

        indffrt=(ns-1)*Pfslen+1:ns*Pfslen;
        resffrt(indffrt,ffrtclnep)=ep; 
        resffrt(indffrt,ffrtclnna)=na; 

        resffrt(indffrt,ffrtclnns)=ns;
        resffrt(indffrt,ffrtclnPsb)=counts(ns).Psb;
        resffrt(indffrt,ffrtclnPsILS)=counts(ns).sucnumILS/Nsamp;
        resffrt(indffrt,ffrtclnPf_req)=counts(ns).Pfs;
        resffrt(indffrt,ffrtclngain)=counts(ns).gain;
        resffrt(indffrt,ffrtclnmu)=counts(ns).muPfs;
        resffrt(indffrt,ffrtclnPfix)=counts(ns).Pfixrts;
        resffrt(indffrt,ffrtclnPscon)=counts(ns).Psconrts;
        
        resffrt(indffrt,ffrtclnPftrue)=counts(ns).Pftruert;
      %  resffrt(indffrt,ffrtclnavails)=counts(ns).muavail; 
        
 
       
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
              Pf2true=counts(ns2).Pftruert(iPfpair);
              Pfix2=counts(ns2).Pfixrts(iPfpair);
              Pscon2=counts(ns2).Psconrts(iPfpair);
              muavails2=counts(ns2).muavail(iPfpair,:);
              
          elseif Pfix1<1
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
                        muavails2=avails2{rows(ns,iPfpair),iPfpair,ns}(muslen,:);
                        Pfix2=0;
                        Pscon2=0;
                    end
            else 
                mu2=1;
                Pf2true=0;
                muavails2=avails2{rows(ns,iPfpair),iPfpair,ns}(muslen,:);
                Pfix2=0;
                Pscon2=0;

            end
            
         availstsrc=counts(ns).muavail(iPfpair,:) + muavails2./Nsamp;
            
         mu1=counts(ns).muPfs(indpf1);
         Pfix1=counts(ns).Pfixrts(indpf1);
         Pscon1=counts(ns).Psconrts(indpf1);
         Pftrue1=counts(ns).Pftruert(indpf1);
         restsrc(indrestsrc,:)=[ep, na, Psb, ns, counts(ns).Psb, counts(ns).gain, Pf1, mu1, Pfix1, Pscon1, Pftrue1, ns2, counts(ns2).Psb, counts(ns2).gain, mu2, Pfix2, Pscon2, Pf2true, availstsrc ];   

%          restsrc=[restsrc; ep, na, Psb, ns, counts(ns).Psb, counts(ns).gain, Pf1, mu1, Pfix1, Pscon1, Pftrue1, ns2, counts(ns2).Psb, counts(ns2).gain, mu2, Pfix2, Pscon2, Pf2true, availstsrc ];   
         %[rowmu21,~,Pf2true]=stdbhat,

        end
    end
