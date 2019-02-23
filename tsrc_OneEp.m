function [restsrc]=tsrc_OneEp(ep,Qa,Nsamp,Psb,PfixData)%Qa,Qb,Qab,Ps, PfixData
%% input: Qa, Qb, Qab, Ps
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
na      = size(Qa,1);
% nb      = size(Qb,1); 
% Qx=[Qa,Qab;
%     Qab', Qb];
% Qx      = (tril(Qx,0)+tril(Qx,-1)');


atrue=randi([-200 200],na,1);% generate random integers
% btrue=zeros(nb,1);
% xtrue=[atrue;btrue];
ncands=2;

[Qzhat,Z,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');


Pfttls=[0.001,0.01]';
Pfttlslen=length(Pfttls);

Pf2s=[0.0001:0.0001:0.001,0.002:0.001:0.01,0.02:0.01:0.1]';
Pf2slen=length(Pf2s);   

Ps0=zeros(Pfttlslen,1);
ns2=Ps0;
for i=1:Pfttlslen
    Pf1=Pfttls(i)*0.9;
    Ps0(i)=1-Pf1;
    ns2(i)=Ps2ns(Psb,Ps0(i),na,D);
end

% exclude the cases when ns2==0.
ns2flag=ns2==0;
ns2=ns2(~ns2flag);
if isempty(ns2), 
    restsrc=zeros(1,12);
    restsrc(1)=ep;
    return;
end
Pfttls=Pfttls(~ns2flag);
Pfttlslen=length(Pfttls);
Pf1s=0.9*Pfttls;

minns2=min(ns2);

mu1s=zeros(Pfttlslen,na);
mu2s=[0.001:0.001:1]';
for ipf=1:Pfttlslen
    Pfreq1=Pf1s(ipf);
    for ns=1:na
        mu1s(ipf,ns)=findmu(ns,Pfreq1,PfixData);% write this function
    end
end
sucflag2=zeros(length(mu2s),Pfttlslen);
rtpass2=zeros(length(mu2s),Pfttlslen);
fixnum2=cell(na,Pfttlslen); 
for i=1:na
    for j=1:Pfttlslen
        fixnum2{i,j}=zeros(length(mu2s),1);
    end
end
%% do the simulation for each sample
sucnum2=fixnum2;
failnum2=fixnum2;
failrate2=fixnum2;
fixrate2=fixnum2;
sucrate2=fixnum2;
unfixnum1=fixnum2;
for i=1:Nsamp
%         xhat=mvnrnd(xtrue,Qx,1)';
%         ahat=xhat(1:na); %bhat=xhat(na+1:na+nb);
%         zhat=Z'*ahat;
        zhat=mvnrnd(ztrue,Qzhat,1)';
        
%         bert=norm(bhat(1:3))*ones(muslen,na);   
%         beflt=norm(bhat(1:3));
%     % Do TSRC here
         for ipf=1:Pfttlslen
             k2=na-ns2(ipf)+1;
             [zpar2,sqnorm2] = ssearch(zhat(k2:end),L(k2:end,k2:end),D(k2:end),ncands);
             sucflag2(:,ipf)=sum(~(ztrue(k2:end)==zpar2(:,1)))==0;
             rtpass2(:,ipf)=sqnorm2(1)/sqnorm2(2)<mu2s;  
         end
         
         for ns=minns2:na
             k1=na-ns+1;
             [zpar1,sqnorm1] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);% I don't care the fix rate and failure rate in the first step currently.
             for ipf=1:Pfttlslen
                 if ns>ns2(ipf)
%                  if ns>=ns2(ipf)
                     if sqnorm1(1)/sqnorm1(2)>mu1s(ipf,ns)% the first step is rejected, entering the second step.
                         
                         unfixnum1{ns,ipf}=unfixnum1{ns,ipf}+1;
                         fixnum2{ns,ipf}=fixnum2{ns,ipf}+rtpass2(:,ipf);
                         sucnum2{ns,ipf}=sucnum2{ns,ipf}+double(sucflag2(:,ipf)&rtpass2(:,ipf));
                         failnum2{ns,ipf}=fixnum2{ns,ipf}-sucnum2{ns,ipf};
                     end
                 end
             end
         end

end

for ns=minns2:na
    for ipf=1:Pfttlslen
        failrate2{ns,ipf}=failnum2{ns,ipf}./unfixnum1{ns,ipf};%fixnum2{ns,ipf};
        fixrate2{ns,ipf}=fixnum2{ns,ipf}./unfixnum1{ns,ipf};%Nsamp;
        sucrate2{ns,ipf}=sucnum2{ns,ipf}./unfixnum1{ns,ipf};%fixnum2{ns,ipf};
        
        nullidx=fixnum2{ns,ipf}==0;
        failrate2{ns,ipf}(nullidx)=1;
        sucrate2{ns,ipf}(nullidx)=0;
        
    end
end
clear fixnum2; clear failnum2; clear sucnum2; clear unfixnum1; 
% calculate the failure rate
themu2=zeros(na,Pf2slen);
thefixrate2=zeros(na,Pf2slen);
thefailrate2=zeros(na,Pf2slen);
thesucrate2=zeros(na,Pf2slen);
ns1slen=0;
for ipf=1:Pfttlslen
    ns1slen=ns1slen+na-ns2(ipf);
end
restsrc=zeros(Pf2slen*ns1slen,12);
resrowcnt=0;
for ipf=1:Pfttlslen
    k2=na-ns2(ipf)+1;
    Pssub2=prod ( 2 * normcdf(1./(2*sqrt(D(k2:end)))) -1 ); 
%     for ns=ns2(ipf):na
    for ns=ns2(ipf)+1:na

        for i=1:Pf2slen
            Pf2req=Pf2s(i);
            [rowmu2,~, Pf2true]=find(failrate2{ns,ipf}<=Pf2req,1,'last');
            if isempty(Pf2true)
                themu2(ns,i)=0;
                thefixrate2(ns,i)=0;
                thefailrate2(ns,i)=0;
                thesucrate2(ns,i)=0;
            else
                themu2(ns,i)=mu2s(rowmu2);
                thefixrate2(ns,i)=fixrate2{ns,ipf}(rowmu2);
                thefailrate2(ns,i)=failrate2{ns,ipf}(rowmu2);
                thesucrate2(ns,i)=sucrate2{ns,ipf}(rowmu2);
            end
        resrowcnt=resrowcnt+1;    
        restsrc(resrowcnt,:)=[ep, na, Psb, Pfttls(ipf), ns, Pf1s(ipf), ns2(ipf), Pssub2, Pf2req, themu2(ns,i),thefixrate2(ns,i), thefailrate2(ns,i)];%Pfix1,

        end
%      indrows=(ipf-1)*(Pf2slen*);
   

    end
end
% display(ep);
% tsrcclnep=1;tsrcclnna=2;tsrcclnPsb=3;tsrcclnPfttl=4;tsrcclnns=5;tsrcclnPf1=6;
% tsrcclnns2=7;tsrcclnPssub2=8;tsrcclnPf2req=9;tsrcclnmu2=10;tsrcclnPfix2=11; tsrcclnPf2=12;
