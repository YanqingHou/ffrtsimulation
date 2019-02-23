function [resffrt,resffdt,restsrc]=muFixRateOneEp(ep,Qa,Qb,Qab,Nsamp,Psb,P0)%Qa,Qb,Qab,Ps
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


counts=struct('ns',[],'Psb',[],'Kcoef1',[],'gain',[],'errs',[],'availsrt',[],'availsdt',[],...
    'ratio',[],'rtmus',[],'diff',[],'dtdiffs',[],'sucnumILS',[],'fltavail',[],...
    'fixnumrt',[],'sucnumrt',[],'failnumrt',[],...
    'fixnumdt',[],'sucnumdt',[],'failnumdt',[],...
    'Pfs',[],'muPfs',[],'DiffPfs',[],...
    'Pfixrts',[],'Psconrts',[],'Pfixdts',[],'Pscondts',[],...
    'muavail',[],'diffavail',[],'Pftruert',[],'Pftruedt',[],...
    'fixnum2',[],'sucnum2',[],'failnum2',[],...
    'mufixnum2',[],'musucnum2',[],'mufailnum2',[]);
% Pss=zeros(n0,1);
mus=0.01:0.01:1;
mus=mus';
muslen=length(mus);
diffs=0:0.1:15;
diffs=diffs';
diffslen=length(diffs);
errs=[0.01:0.01:0.1,0.2:0.1:0.5]';
errslen=length(errs);
Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';
Pfslen=length(Pfs);
% avails=zeros(muslen,errslen);

if n0==0 
    resffrt=zeros(1,11+errslen);
    resffrt(1)=ep;
    
    resffdt=zeros(1,11+errslen);
    resffrt(1)=ep;
    
    restsrc=zeros(1,19);
    restsrc(1)=ep;    
    return ; 
end
resffrt=zeros(n0*(1+Pfslen),11+errslen);% one more row for float solution
resffdt=zeros(n0*(1+Pfslen),11+errslen);
restsrc=zeros(n0*Pfslen*muslen,19);

avails2=cell(muslen,Pfslen,n0);
Ps0=zeros(1,Pfslen);
ns2=Ps0;
for i=1:Pfslen
Ps0(i)=1-Pfs(i);
ns2(i)=Ps2ns(Psb,Ps0(i),na,D);
end
% initialization of the struct array
for ns=1:n0
    k1=na-ns+1;
%     Pss(ns)= prod (2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
    counts(ns).ns=ns;
    counts(ns).Psb=prod(2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
    counts(ns).sucnumILS=0;
    counts(ns).rtmus=mus;
    counts(ns).fltavail = zeros(1,errslen); 
    
    Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
    Kcoef1=Qbzpar1/Qzpar1;  
    Qbp1=Qb-Kcoef1*Qbzpar1';
    Rs1=sqrt(det(Qb(1:3,1:3))/det(Qbp1(1:3,1:3)))^(1/3);%large than 1
    counts(ns).gain=Rs1;
    counts(ns).Kcoef1=Kcoef1;
    
    counts(ns).fixnumrt=zeros(muslen,1);
    counts(ns).sucnumrt=counts(ns).fixnumrt;
    counts(ns).failnumrt=counts(ns).fixnumrt;
    
    counts(ns).dtdiffs=diffs;
    counts(ns).fixnumdt=zeros(diffslen,1);
    counts(ns).sucnumdt=counts(ns).fixnumdt;
    counts(ns).failnumdt=counts(ns).fixnumdt;
    

    counts(ns).errs=errs;  % for saving in the output file
    counts(ns).availsrt=zeros(muslen,errslen); % for each mu and each required err, there is an availability
    counts(ns).availsdt=zeros(diffslen,errslen);
    counts(ns).muavail=zeros(Pfslen,errslen);  % There are 11 Pfs, for each required Pf, find the mu value; and for this mu value, calculate the availability of 20 errors.
    counts(ns).diffavail=zeros(Pfslen,errslen);

    
    counts(ns).Pfs=Pfs;
    counts(ns).fixnum2=zeros(muslen,muslen,Pfslen);
    counts(ns).sucnum2=zeros(muslen,muslen,Pfslen);
    counts(ns).failnum2=zeros(muslen,muslen,Pfslen);
    
    counts(ns).musucnum2=zeros(Pfslen,muslen);
    counts(ns).mufixnum2=zeros(Pfslen,muslen);
    counts(ns).mufailnum2=zeros(Pfslen,muslen);

    counts(ns).muPfs=zeros(Pfslen,1);
    counts(ns).DiffPfs=zeros(Pfslen,1);
    counts(ns).Pfixrts=zeros(Pfslen,1);
    counts(ns).Psconrts=zeros(Pfslen,1);
    counts(ns).Pfixdts=zeros(Pfslen,1);

    counts(ns).Pscondts=zeros(Pfslen,1);
    counts(ns).Pftruert=zeros(Pfslen,1);
    counts(ns).Pftruedt=zeros(Pfslen,1);
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
        bedt=bert;
        
%         berttsrc=bert;   

%         zhat=mvnrnd(ztrue,Qzhat,1)';
        sucflags=zeros(n0,1);
    for ns=1:n0
        k1=na-ns+1; %zp=zhat(na-ns+1:na);
        [zpar,sqnorm] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);
        zpar=zpar(:,1);
        bchk=bhat-counts(ns).Kcoef1*(zhat(k1:end)-zpar);
        
        sucflag=sum(~(ztrue(k1:end)==zpar))==0;
        counts(ns).sucnumILS = counts(ns).sucnumILS+sucflag;         % Correctly fixed
        sucflags(ns)=sucflag;
        counts(ns).fltavail = counts(ns).fltavail+double(bert(1)<=counts(ns).errs');
        counts(ns).ratio=sqnorm(1)/sqnorm(2);
        
        rtpass=sqnorm(1)/sqnorm(2)<counts(ns).rtmus;                % 100 different mu
        counts(ns).fixnumrt=counts(ns).fixnumrt+rtpass;             % accepted by FFRT
        counts(ns).sucnumrt=counts(ns).sucnumrt + double(rtpass&sucflag);   % success rate after FFRT: (Accepted by FFRT & correct) / N
        counts(ns).failnumrt=counts(ns).fixnumrt-counts(ns).sucnumrt;% failure rate after FFRT (Fix-Fix & correct)
        bert(rtpass,ns)=norm(bchk(1:3));                                % if FFRT passed, baseline estimation is updated, and the rms is updated.
        resmat=comparevect(bert(:,ns),counts(ns).errs);                    % compare the rms with required errors for each mu in FFRT, and this result in a matrix.
        counts(ns).availsrt=counts(ns).availsrt+resmat;                 % In the end, it needs to be devided by N---- the total sample number.
        counts(ns).diff=sqnorm(2)-sqnorm(1);
        
        dtpass=sqnorm(2)-sqnorm(1)>=counts(ns).dtdiffs;
        counts(ns).fixnumdt=counts(ns).fixnumdt+dtpass;             % accepted by FFDT
        counts(ns).sucnumdt=counts(ns).sucnumdt + double(dtpass&sucflag);   % success rate after FFDT: (Accepted by FFDT & correct) / N
        counts(ns).failnumdt=counts(ns).fixnumdt-counts(ns).sucnumdt;% failure rate after FFDT (Fix-Fix & correct)
                
        bedt(dtpass,ns)=norm(bchk(1:3));                                % if FFRT passed, baseline estimation is updated, and the rms is updated.
        resmat=comparevect(bedt(:,ns),counts(ns).errs);                    
        counts(ns).availsdt=counts(ns).availsdt+resmat;           
    end
    % Do TSRC here
   for j=1:Pfslen % for every Pf... wait, is it neccessary?
       rtpass2=counts(ns2(j)).ratio<=counts(ns2(j)).rtmus; % the ratio test of the second step
        for ns=ns2(j)+1:n0
            rtpass=counts(ns).ratio<=counts(ns).rtmus;                % 100 different mu
            rtreject=~rtpass';
%             for k=1:length(rtpass)
%                 if rtpass(k)==0                                      % the first step of TSRC rejected, go to the second step                   
                    counts(ns).fixnum2(:,:,j)=counts(ns).fixnum2(:,:,j)+kron(rtpass2,rtreject); % the first index is for mus in step2, the second index is for mus in first step, the third index is for different Pf                                       %
                    counts(ns).sucnum2(:,:,j)=counts(ns).sucnum2(:,:,j) + double( kron(rtpass2,rtreject)&sucflags(ns2(j)));   % success rate after FFRT: (Accepted by FFRT & correct) / N
                    counts(ns).failnum2(:,:,j)=counts(ns).failnum2(:,:,j)-counts(ns).sucnum2(:,:,j); 
                    for k=1:muslen  % for each mu in the first step of TSRC
                        if rtreject(k)==1
                           resmat=comparevect(bert(:,ns2(j)),counts(ns).errs);                    % compare the rms with required errors for each mu in FFRT, and this result in a matrix.
                           avails2{k,j,ns}=avails2{k,j,ns}+resmat;% float 和 fixed solution 要分开加
                        end
                    end
%                     avails2=cell(muslen,Pfslen,n0);
%                     avails2{j,i,ns}=zeros(muslen,errslen);
% 
%                     counts(ns).avails2(:,:,j)=counts(ns).avails2(:,:,j)
%                 end
%             end

        end
   end
end

%all the possible result simulated. In the next, we need to find the best
%parameter for FFRT and TSRC

ffrtclnep=1;ffrtclnna=2;ffrtclnns=3;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;
ffrtclngain=7;ffrtclnmu=8;ffrtclnPfix=9;ffrtclnPscon=10;
ffrtclnPftrue=11;ffrtclnavails=12:12+errslen-1;

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
         0.01, 0.009];

tsrcclnep=1; tsrcclnna=2; tsrcclnPsb=3; tsrcclnstdbhat=4;
tsrcclnns1=5; tsrcclnPssub1=6; tsrcclngain1=7; tsrcclnPf_req1=8; tsrcclnmu1=9; tsrcclnPfix1=10;tsrcclnPscon1=11; tsrcclnPftrue1=12;
tsrcclnns2=13; tsrcclnPssub2=14; tsrcclngain2=15; tsrcclnmu2=16; tsrcclnPfix2=17;tsrcclnPscon2=18; tsrcclnPftrue2=19; 
tsrcclnavails=20:20+errslen-1;

% parameter for FFRT
for ns=1:n0
    
    counts(ns).sucnumrt=counts(ns).sucnumrt./counts(ns).fixnumrt;
    counts(ns).failnumrt=counts(ns).failnumrt/Nsamp;
    counts(ns).fixnumrt=counts(ns).fixnumrt/Nsamp;
    counts(ns).availsrt=counts(ns).availsrt/Nsamp;
    
    counts(ns).sucnumdt=counts(ns).sucnumdt./counts(ns).fixnumdt;
    counts(ns).failnumdt=counts(ns).failnumdt/Nsamp;
    counts(ns).fixnumdt=counts(ns).fixnumdt/Nsamp;
    counts(ns).availsdt=counts(ns).availsdt/Nsamp;
    
    row=zeros(Pfslen,1);
    rowdt=row;
    for pfi=1:Pfslen
        [row(pfi), ~, mupf]=find(counts(ns).failnumrt<=Pfs(pfi),1,'last'); % find the critical value with a given failure rate
        [rowdt(pfi),~, DiffPf1]=find(counts(ns).failnumdt<=Pfs(pfi),1,'first');
        if isempty(mupf)
%            mupf=0; 
           counts(ns).muPfs(pfi)=0;%mupf;
           counts(ns).Pfixrts(pfi)=0;%counts(ns).fixnumrt(row);
           counts(ns).Psconrts(pfi)=0;%counts(ns).sucnumrt(row);
           counts(ns).muavail(pfi,:)=counts(ns).availsrt(muslen,:);%
           counts(ns).Pftruert(pfi)=0;%counts(ns).failnumrt(row);
       
           counts(ns).musucnum2(pfi,:)=zeros(1,muslen);%for TSRC
           counts(ns).mufixnum2(pfi,:)=zeros(1,muslen);%for TSRC
           counts(ns).mufailnum2(pfi,:)=zeros(1,muslen);%for TSRC
        else
            counts(ns).muPfs(pfi)=mupf;
            counts(ns).Pfixrts(pfi)=counts(ns).fixnumrt(row(pfi));
            counts(ns).Psconrts(pfi)=counts(ns).sucnumrt(row(pfi));
            counts(ns).muavail(pfi,:)=counts(ns).availsrt(row(pfi),:);
            counts(ns).Pftruert(pfi)=counts(ns).failnumrt(row(pfi));
            
            %for TSRC
            counts(ns).musucnum2(pfi,:)=counts(ns).sucnum2(:,row(pfi),pfi)./counts(ns).fixnum2(:,row(pfi),pfi); nanidx=isnan(counts(ns).musucnum2(pfi,:)); counts(ns).musucnum2(pfi,nanidx)=0; 
            counts(ns).mufixnum2(pfi,:)=counts(ns).fixnum2(:,row(pfi),pfi)./(Nsamp*(1-counts(ns).Pfixrts(pfi))); nanidx=isnan(counts(ns).mufixnum2(pfi,:)); counts(ns).mufixnum2(pfi,nanidx)=0; 
            counts(ns).mufailnum2(pfi,:)=counts(ns).failnum2(:,row(pfi),pfi)./counts(ns).fixnum2(:,row(pfi),pfi); nanidx=isnan(counts(ns).mufailnum2(pfi,:)); counts(ns).mufailnum2(pfi,nanidx)=1;
             % 在这种情况下，成功率和固定率可以赋值为0，但是失败率赋值为0就不对了。
%             % When Pfix is nonzero, or when Pf2req 
%             Pf2req=(Pfttl-counts(ns).Pftruert(pfi))/(1-counts(ns).Pfixrts(pfi));
%             mu2=find(counts(ns).mufailnum2(pfi,:)<=Pf2req,1,'last');
            
        end
         
          if isempty(DiffPf1)
%            mupf=0; 
           counts(ns).DiffPfs(pfi)=1000;%mupf;
           counts(ns).Pfixdts(pfi)=0;%counts(ns).fixnumrt(row);
           counts(ns).Pscondts(pfi)=0;%counts(ns).sucnumrt(row);
           counts(ns).diffavail(pfi,:)=counts(ns).availsdt(muslen,:);%
           counts(ns).Pftruedt(pfi)=0;%counts(ns).failnumrt(row);
      
        else
            counts(ns).DiffPfs(pfi)=DiffPf1;
            counts(ns).Pfixdts(pfi)=counts(ns).fixnumdt(rowdt(pfi));
            counts(ns).Pscondts(pfi)=counts(ns).sucnumdt(rowdt(pfi));
            counts(ns).diffavail(pfi,:)=counts(ns).availsdt(rowdt(pfi),:);
            counts(ns).Pftruedt(pfi)=counts(ns).failnumdt(rowdt(pfi));
         end
    end

%     resffrt(ns,:)=[ep,counts(ns).Pftruert(pfi),];

        indffrt=(ns-1)*(Pfslen+1)+2:ns*(Pfslen+1);
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
        resffrt(indffrt,ffrtclnavails)=counts(ns).muavail; 
        
        
        indflt=(ns-1)*(Pfslen+1)+1; % for float solution
        
        resffrt(indffrt,ffrtclnep)=ep; 
        resffrt(indffrt,ffrtclnna)=na; 

        resffrt(indffrt,ffrtclnns)=ns;
        resffrt(indffrt,ffrtclnPsb)=counts(ns).Psb;
        resffrt(indffrt,ffrtclnPsILS)=counts(ns).sucnumILS/Nsamp;
        
        resffrt(indflt,ffrtclnPf_req)=1;
        resffrt(indflt,ffrtclngain)=1;
        resffrt(indflt,ffrtclnmu)=1;
        resffrt(indflt,ffrtclnPfix)=0;
        resffrt(indffrt,ffrtclnPscon)=0;
        resffrt(indflt,ffrtclnPftrue)=1;
        resffrt(indflt,ffrtclnavails)=counts(ns).fltavail;
        
        
        indtsrc=(ns-1)*(Pfslen)*muslen+1:ns*(Pfslen)*muslen;
        restsrc(indtsrc,tsrcclnep)=ep;
        restsrc(indtsrc,tsrcclnna)=na;
        restsrc(indtsrc,tsrcclnns1)=ns;
        restsrc(indtsrc,tsrcclnPsb)=Psb;
        restsrc(indtsrc,tsrcclnPssub1)=counts(ns).Psb;
        restsrc(indtsrc,tsrcclngain1)=counts(ns).gain;
        for j=1:Pfslen
            indtsrc2=(ns-1)*(Pfslen)*muslen+(j-1)*muslen+1:(ns-1)*(Pfslen)*muslen+j*muslen;
            restsrc(indtsrc2,tsrcclnPf_req1)=counts(ns).Pfs(j);
            restsrc(indtsrc2,tsrcclnmu1)=counts(ns).muPfs(j);
            restsrc(indtsrc2,tsrcclnPfix1)=counts(ns).Pfixrts(j);
            restsrc(indtsrc2,tsrcclnPscon1)=counts(ns).Psconrts(j);
            restsrc(indtsrc2,tsrcclnPftrue1)=counts(ns).Pftruert(j);
            
            restsrc(indtsrc2,tsrcclnns2)=ns2(j);
            restsrc(indtsrc2,tsrcclnPssub2)=counts(ns2(j)).Psb;
            restsrc(indtsrc2,tsrcclngain2)=counts(ns2(j)).gain;
            restsrc(indtsrc2,tsrcclnmu2)=mus;
            restsrc(indtsrc2,tsrcclnPfix2)=counts(ns).mufixnum2(j,:)';
            restsrc(indtsrc2,tsrcclnPscon2)=counts(ns).musucnum2(j,:)';
            restsrc(indtsrc2,tsrcclnPftrue2)=counts(ns).mufailnum2(j,:)';
%             Pf2req=(Pfreq-Pfreq1)/(1-Pfix1);
% avails{k,j,ns}counts(ns).muavail(:,:,j)+
            restsrc(indtsrc2,tsrcclnavails)=kron(ones(muslen,1),counts(ns).muavail(j,:))+avails2{row(j),j,ns};%counts(ns2(j)).muavail(j,:);

         end

        
        indffdt=(ns-1)*(Pfslen+1)+2:ns*(Pfslen+1);
        resffdt(indffdt,ffrtclnep)=ep; 
        resffdt(indffdt,ffrtclnna)=na;

        resffdt(indffdt,ffrtclnns)=ns;
        resffdt(indffdt,ffrtclnPsb)=counts(ns).Psb;
        resffdt(indffdt,ffrtclnPsILS)=counts(ns).sucnumILS/Nsamp;
        resffdt(indffdt,ffrtclnPf_req)=counts(ns).Pfs;
        resffdt(indffdt,ffrtclngain)=counts(ns).gain;
        resffdt(indffdt,ffrtclndiff)=counts(ns).DiffPfs;
        resffdt(indffdt,ffrtclnPfix)=counts(ns).Pfixdts;
        resffdt(indffdt,ffrtclnPscon)=counts(ns).Pscondts;

        resffdt(indffdt,ffrtclnPftrue)=counts(ns).Pftruedt;
        resffdt(indffdt,ffrtclnavails)=counts(ns).diffavail;
        
        indflt=(ns-1)*(Pfslen+1)+1; % for float solution
        
        resffdt(indffdt,ffrtclnep)=ep; 
        resffdt(indffdt,ffrtclnna)=na;

        resffdt(indffdt,ffrtclnns)=ns;
        resffdt(indffdt,ffrtclnPsb)=counts(ns).Psb;
        resffdt(indffdt,ffrtclnPsILS)=counts(ns).sucnumILS/Nsamp;
        resffdt(indflt,ffrtclnPf_req)=1;
        resffdt(indflt,ffrtclngain)=1;
        resffdt(indflt,ffrtclnmu)=1;
        resffdt(indflt,ffrtclnPfix)=0;
        resffdt(indffdt,ffrtclnPscon)=0;

        resffdt(indflt,ffrtclnPftrue)=1;
        resffdt(indflt,ffrtclnavails)=counts(ns).fltavail;
end

    %for TSRC only.等下要把上一个循环中关于TSRC的部分删除掉，只留FFRT的部分。
for ns=ns2+1:na
    for iPfpair=1:size(Pfpairs,1)
        Pfttl=Pfpairs(iPfpair,1);
        Pf1=Pfpairs(iPfpair,2);
        
         restsrc(indtsrc2,tsrcclnPf_req1)=Pf1;%counts(ns).Pfs(j);
         restsrc(indtsrc2,tsrcclnmu1)=counts(ns).muPfs(iPfpair);
         Pfix1=counts(ns).Pfixrts(iPfpair);
         Pf2=(1-Pf1)/(1-Pfix1);
         [iPfpair2,Pf21]=find(Pfpairs(:,2)<=Pf2,1,'last');
%                      restsrc(indtsrc2,tsrcclnPftrue2)=
         [rowmu2, mu2]=find(counts(ns).mufailnum2(iPfpair2,:)'<=Pf2,1,'last');
         Pfix2=counts(ns).mufixnum2(iPfpair2,rowmu2);
         Pf2true=counts(ns).mufailnum2(iPfpair2,rowmu2);
         restsrc(indtsrc2,tsrcclnPfix1)=counts(ns).Pfixrts(j);
         restsrc(indtsrc2,tsrcclnPscon1)=counts(ns).Psconrts(j);
         restsrc(indtsrc2,tsrcclnPftrue1)=counts(ns).Pftruert(j);
            
        %[ep, na, Psb, ns1, Pssub1, gain1, Pf_ttl_req, Pf1_req, mu1, Pfix1, Pscon1, Pftrue1, ns2, Pssub2, gain2, Pf_req2, mu2, Pfix2, Pftrue2]
    end
end


