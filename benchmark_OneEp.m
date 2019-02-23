function [resbenchmark]=benchmark_OneEp(ep,Qa,Qb,Qab,Nsamp,Psb,PfixData,TSRCData)
%% Input: 1. PfixData
%         2. TSRCData
%% 0.Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
na      = size(Qa,1);
nb      = size(Qb,1);
Qx=[Qa,Qab;
    Qab', Qb];
Qx      = (tril(Qx,0)+tril(Qx,-1)');


atrue=randi([-200 200],na,1);% generate random integers
btrue=zeros(nb,1);
xtrue=[atrue;btrue];
ncands=2;

[Qzhat,Z,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');
%% 1.the above is initialized infomation, not need to double check%, this section do more initialization

Pfttls=[0.001,0.01]';
Pfttlslen=length(Pfttls);

Ps0=zeros(Pfttlslen,1);
ns2=Ps0;
g2=Ps0;
for i=1:Pfttlslen
    Pf1=Pfttls(i)*0.9;
    Ps0(i)=1-Pf1;
    ns2(i)=Ps2ns(Psb,Ps0(i),na,D);
    
    k2=na-ns2(i)+1;
    Qzpar2 = Qzhat(k2:end,k2:end); Zpar2 = Z(:,k2:end); Qbzpar2=Qab'*Zpar2;
    Kcoef2=Qbzpar2/Qzpar2;
    Qbp2=Qb-Kcoef2*Qbzpar2';
    g2(i)=sqrt(det(Qb(1:3,1:3))/det(Qbp2(1:3,1:3)))^(1/3);%large than 1
    
end

Pf1s=0.9*Pfttls;

mu1s=zeros(Pfttlslen,1);
optsub=zeros(Pfttlslen,1);
mu2s=zeros(Pfttlslen,1);
Kcoef1s=cell(Pfttlslen,1);
Kcoef2s=cell(Pfttlslen,1);
Pfix1s=zeros(Pfttlslen,1);
Pf1trues=ones(Pfttlslen,1);
%
% optsubfit=optsub;
% Kcoef1fits=Kcoef1s;
% Pfix1fits=Pfix1s;
% mufit1s=Pfix1s;
% tmp=load('MuFit_exp_coef.mat');
% fitresmus=tmp.fitresmus; clear tmp;
% tmp=load('PfixFit_exp_coef.mat');
% fitresPfixs=tmp.fitresPfixs; clear tmp;
for ipf=1:Pfttlslen
    Pfreq1=Pf1s(ipf);
    %     find the optimal subset that maximizes Pfix*(g1-g2)
    maxdEg=0;
    %     maxdEgfit=0;
    optsub(ipf)=ns2(ipf);
    if ns2(ipf)>0
        k1=na-ns2(ipf)+1;
        Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
        Kcoef1=Qbzpar1/Qzpar1;
        Kcoef1s{ipf,1}=Kcoef1;
    end
    for ns=ns2(ipf)+1:na
        k1=na-ns+1;
        %         Ps1= prod (2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
        %         Pfix1=cpfixrate(Ps1,Pf,ns);
        
        Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
        Kcoef1=Qbzpar1/Qzpar1;
        Qbp1=Qb-Kcoef1*Qbzpar1';
        g1=sqrt(det(Qb(1:3,1:3))/det(Qbp1(1:3,1:3)))^(1/3);%large than 1
        
        Pfix1=findPfix(ns,Pfreq1,PfixData);%%%%%%%%%%%%%%write this function
        %         cfunc=fitresPfixs{ns,ipf};
        %         Psbns= prod ( 2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
        %         Pfix1fit1=feval(cfunc,1-Psbns);
        %         cfunc=fitresmus{ns,ipf};
        %         mu1fit=feval(cfunc,1-Psbns);
        %         Pfix1fit2=
        %         Pfix1fit3=
        dEg=Pfix1*(g1-g2(ipf));
        
        %         dEgfit=Pfix1fit*(g1-g2(ipf));
        %         if dEgfit>maxdEgfit
        %             maxdEgfit=dEgfit;
        %             optsubfit(ipf)=ns;
        %             Kcoef1fits{ipf,1}=Kcoef1;
        %             Pfix1fits(ipf)=Pfix1fit1;
        %             mufit1s(ipf)=mu1fit;
        %         end
        
        
        if dEg>maxdEg
            maxdEg=dEg;
            optsub(ipf)=ns;
            Kcoef1s{ipf,1}=Kcoef1;
            Pfix1s(ipf)=Pfix1;
            %             Pf1trues(ipf)=Pf1true;
        end
    end
    %     find optimal subset
    %     find the mu1
    [mu1s(ipf), Pf1trues(ipf)]=findmu(optsub(ipf),Pfreq1,PfixData);% function written
    
    if ns2(ipf)>0
        if Pfix1s(ipf)==0
            Pf2req=Pfttls(ipf);
        else
%            Pf2req=(Pfttls(ipf)-Pf1s(ipf))/(1-Pfix1s(ipf));
            Pf2req=(Pfttls(ipf)-Pf1trues(ipf))/(1-Pfix1s(ipf));
            Ns2=Nsamp*(1-Pfix1s(ipf));
            deci=100/Ns2;
  
            tmpdec1=mod(Pf2req,deci);
            Pf2req=Pf2req-tmpdec1;
           % if Pf2req<deci
            %    Pf2req=0;
           % end
         end
        %           mu2=findmu2(Pfttl,ns1, ns2,Pf2req, TSRCData)
        mu2s(ipf)=findmu2(Pfttls(ipf),optsub(ipf),ns2(ipf),Pf2req,TSRCData);% write this function
        k2=na-ns2(ipf)+1;
        Qzpar2 = Qzhat(k2:end,k2:end); Zpar2 = Z(:,k2:end); Qbzpar2=Qab'*Zpar2;
        Kcoef2=Qbzpar2/Qzpar2;
        Kcoef2s{ipf,1}=Kcoef2;
    end
end

fixnum2=zeros(Pfttlslen,1);
sucnum2=fixnum2;
failnum2=fixnum2;
fixnum1=fixnum2;
sucnum1=fixnum2;
failnum1=fixnum2;

errs=[0.01:0.01:0.1,0.2:0.1:0.5];
errslen=length(errs);
availnum=zeros(Pfttlslen,errslen);
% availnumfit=availnum;
%% 2.The calculation of N samples.
for i=1:Nsamp
    xhat=mvnrnd(xtrue,Qx,1)';
    ahat=xhat(1:na); bhat=xhat(na+1:na+nb);
    zhat=Z'*ahat;
    beflt=norm(bhat(1:3));
    
    for ipf=1:Pfttlslen
        if optsub(ipf)>0
            k1=na-optsub(ipf)+1;
            [zpar1,sqnorm1] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);
            sucflag1=sum(~(ztrue(k1:end)==zpar1(:,1)))==0;
            if sqnorm1(1)/sqnorm1(2)>mu1s(ipf)
                if ns2(ipf)>0
                    k2=na-ns2(ipf)+1;
                    [zpar2,sqnorm2] = ssearch(zhat(k2:end),L(k2:end,k2:end),D(k2:end),ncands);
                    sucflag2=sum(~(ztrue(k2:end)==zpar2(:,1)))==0;
                    %                     ratio test for the second step.
                    if sqnorm2(1)/sqnorm2(2)>mu2s(ipf)
                        %                             use the float solution, update the
                        %                             availability
                        availnum(ipf,:)=availnum(ipf,:)+double(beflt<=errs);
                        
                    else
                        fixnum2(ipf)=fixnum2(ipf)+1;
                        sucnum2(ipf)=sucnum2(ipf)+sucflag2;
                        failnum2(ipf)=fixnum2(ipf)-sucnum2(ipf);
                        bchk2=bhat-Kcoef2s{ipf,1}*(zhat(k2:end)-zpar2(:,1));
                        availnum(ipf,:)=availnum(ipf,:)+double(norm(bchk2(1:3))<=errs);
                        %                             update the baseline solution with ns2(ipf)
                    end
                else % ns2==0, use float solution.
                    availnum(ipf,:)=availnum(ipf,:)+double(beflt<=errs);
                end
            else % subset in the first step accepted
                fixnum1(ipf)=fixnum1(ipf)+1;
                sucnum1(ipf)=sucnum1(ipf)+sucflag1;
                failnum1(ipf)=fixnum1(ipf)-sucnum1(ipf);
                bchk1=bhat-Kcoef1s{ipf,1}*(zhat(k1:end)-zpar1(:,1));
                availnum(ipf,:)=availnum(ipf,:)+double(norm(bchk1(1:3))<=errs);
                %                       use the optimal subset
                %                       collect the fix num, success num and fail num
                %                       calculate the availability
            end
            
        else % optsub(ipf)==0
            %                   use the float solution
            availnum(ipf,:)=availnum(ipf,:)+double(beflt<=errs);
        end
    end
end

fixrate1=fixnum1/Nsamp;
sucrate1=sucnum1/Nsamp;
failrate1=failnum1/Nsamp;
fixrate2=fixnum2./(Nsamp-fixnum1);
sucrate2=sucnum2./(Nsamp-fixnum1);
failrate2=failnum2./(Nsamp-fixnum1);

fixratettl=(fixnum1+fixnum2)/Nsamp;
sucratettl=(sucnum1+sucnum2)/Nsamp;
failratettl=(failnum1+failnum2)/Nsamp;
availability=availnum/Nsamp;
% resbenchmark should also contain required total failure rate.
resbenchmark=[ep*ones(Pfttlslen,1), Pfttls, fixrate1, sucrate1, failrate1,fixrate2,sucrate2,failrate2, fixratettl, sucratettl, failratettl, availability];
nanidx=isnan(resbenchmark);
resbenchmark(nanidx)=0;
display(ep);
%%End of the file-0222220**********************9
