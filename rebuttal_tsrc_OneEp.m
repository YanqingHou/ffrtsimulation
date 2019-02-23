function [resallARmethods]=rebuttal_tsrc_OneEp(ep,Qa,Qb,Qab,Nsamp,Psb,PfixData,TSRCData)
%% Input: 1. PfixData
%         2. TSRCData
%         3. F_mu, F_fix1, F_mu2
%% 0.Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
na      = size(Qa,1);
nb      = size(Qb,1);
Qx=[Qa,Qab;
    Qab', Qb];
Qx      = (tril(Qx,0)+tril(Qx,-1)');

atrue=randi([-200 200],na,1);% generate random integers
btrue=zeros(nb,1);
xtrue=[atrue;btrue];

[Qzhat,Z,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');

Qbz=Qab'*Z;
Kfar=Qbz/Qzhat;



%% 1.the above is initialized infomation, not need to double check%, this section do more initialization

Pfttls=[0.001,0.01]';
Pfttlslen=length(Pfttls);

Ps0=zeros(Pfttlslen,1);
ns2=Ps0;
g2=Ps0;

% mu1s=zeros(Pfttlslen,1);      Fitmu1s=mu1s;            FARmu1s=mu1s;
optsub=zeros(Pfttlslen,1);    Fitoptsub=optsub;
% mu2s=zeros(Pfttlslen,1);      Fitmu2s=mu2s;
% Kcoef1s=cell(Pfttlslen,1);    FitKcoef1s=Kcoef1s;      FARKcoef=0;
% Kcoef2s=cell(Pfttlslen,1);
% Pfix1s=zeros(Pfttlslen,1);    FitPfix1s=Pfix1s;
% Pf1trues=ones(Pfttlslen,1);   FitPf1trues=Pf1trues;    FARPftrue=Pf1trues;    SRCPftrues=Pf1trues;
Pf1s=Pfttls;%0.9*
FARmus=zeros(Pfttlslen,1);

for i=1:Pfttlslen
    Pf1=Pfttls(i);%*0.9;
    Ps0(i)=1-Pf1;
    ns2(i)=Ps2ns(Psb,Ps0(i),na,D);
    
    k2=na-ns2(i)+1;
    Qzpar2 = Qzhat(k2:end,k2:end); Zpar2 = Z(:,k2:end); Qbzpar2=Qab'*Zpar2;
    Kcoef2=Qbzpar2/Qzpar2;
    Qbp2=Qb-Kcoef2*Qbzpar2';
    g2(i)=sqrt(det(Qb(1:3,1:3))/det(Qbp2(1:3,1:3)))^(1/3);%large than 1
    
    FARmus(i)=cpmu1(1-Psb,Pfttls(i),na);
end


DoStep2flags=zeros(Pfttlslen,1); FitDoStep2flags=DoStep2flags;
paras=struct('mu1',[],'mu2',[],'K1',[],'K2',[]); Fitparas=paras;
for ipf=1:Pfttlslen
    Pfreq1=Pf1s(ipf);
    %     find the optimal subset that maximizes Pfix*(g1-g2)
    EstPfixMethod=0;% not use fitting function
    [optsub(ipf),DoStep2flags(ipf),paras(ipf)]=SearchOptsub(Pfttls(ipf),Pfreq1,Qzhat,Z,D,Qab,Qb,ns2(ipf),g2(ipf),PfixData,TSRCData,EstPfixMethod);
   
    EstPfixMethod=1;% use fitting function
    [Fitoptsub(ipf),FitDoStep2flags(ipf),Fitparas(ipf)]=SearchOptsub(Pfttls(ipf),Pfreq1,Qzhat,Z,D,Qab,Qb,ns2(ipf),g2(ipf),PfixData,TSRCData,EstPfixMethod);
   
end

fixnum2=zeros(Pfttlslen,1);  Fitfixnum2=fixnum2;  SRCfixnum2=fixnum2;  FARfixnum2=fixnum2;
sucnum2=fixnum2;             Fitsucnum2=fixnum2;  SRCsucnum2=sucnum2;  FARsucnum2=sucnum2;
% failnum2=fixnum2;            Fitfailnum2=fixnum2;
fixnum1=fixnum2;             Fitfixnum1=fixnum2;  SRCfixnum1=fixnum1;  FARfixnum1=fixnum1;
sucnum1=fixnum2;             Fitsucnum1=fixnum2;  SRCsucnum1=sucnum1;  FARsucnum1=sucnum1;
% failnum1=fixnum2;            Fitfailnum1=fixnum2;

errs=[0.01:0.01:0.1,0.2:0.1:0.5];
errslen=length(errs);
availnum=zeros(Pfttlslen,errslen);  Fitavailnum=availnum; SRCavailnum=availnum; FARavailnum=availnum;

% availnumfit=availnum;
%% 2.The calculation of N samples.
for i=1:Nsamp
    xhat=mvnrnd(xtrue,Qx,1)';
    ahat=xhat(1:na); bhat=xhat(na+1:na+nb);
    zhat=Z'*ahat;
    beflt=norm(bhat(1:3));
    
    for ipf=1:Pfttlslen
%           [FARavailnum(ipf,:),FARfixnum2(ipf),FARsucnum2(ipf),FARfixnum1(ipf),FARsucnum1(ipf)]=DoFAR(zhat,ztrue,bhat,L,D,...
%                                                                              FARmus(ipf),Kfar,beflt,errs,...
%                                                                              FARavailnum(ipf,:),FARfixnum2(ipf),FARsucnum2(ipf),FARfixnum1(ipf),FARsucnum1(ipf));    
%         if ns2(ipf)==na%
%             
%             availnum(ipf,:)=FARavailnum(ipf,:); fixnum2(ipf)=FARfixnum2(ipf);sucnum2(ipf)=FARsucnum2(ipf);fixnum1(ipf)=FARfixnum1(ipf);sucnum1(ipf)=FARsucnum1(ipf);
%             Fitavailnum(ipf,:)=FARavailnum(ipf,:); Fitfixnum2(ipf)=FARfixnum2(ipf);Fitsucnum2(ipf)=FARsucnum2(ipf);Fitfixnum1(ipf)=FARfixnum1(ipf);Fitsucnum1(ipf)=FARsucnum1(ipf);
%             SRCavailnum(ipf,:)=FARavailnum(ipf,:); SRCfixnum2(ipf)=FARfixnum1(ipf);SRCsucnum2(ipf)=FARsucnum1(ipf);SRCfixnum1(ipf)=FARfixnum2(ipf);SRCsucnum1(ipf)=FARsucnum2(ipf);
%         else
        [availnum(ipf,:),fixnum2(ipf),sucnum2(ipf),fixnum1(ipf),sucnum1(ipf)]=DoTSRCstep1(zhat,ztrue,bhat,L,D,na,optsub(ipf),ns2(ipf),paras(ipf).mu1,paras(ipf).mu2,...
            paras(ipf).K1,paras(ipf).K2,beflt,errs,...
            availnum(ipf,:),fixnum2(ipf),sucnum2(ipf),fixnum1(ipf),sucnum1(ipf));
       
       [Fitavailnum(ipf,:),Fitfixnum2(ipf),Fitsucnum2(ipf),Fitfixnum1(ipf),Fitsucnum1(ipf)]=DoTSRCstep1(zhat,ztrue,bhat,L,D,na,Fitoptsub(ipf),ns2(ipf),Fitparas(ipf).mu1,Fitparas(ipf).mu2,...
            Fitparas(ipf).K1,Fitparas(ipf).K2,beflt,errs,...
            Fitavailnum(ipf,:),Fitfixnum2(ipf),Fitsucnum2(ipf),Fitfixnum1(ipf),Fitsucnum1(ipf));
       
%         [SRCavailnum(ipf,:),SRCfixnum2(ipf),SRCsucnum2(ipf),SRCfixnum1(ipf),SRCsucnum1(ipf)]=DoSRC(zhat,ztrue,bhat,L,D,na,ns2(ipf),...
%                                                                              Fitparas(ipf).K2,beflt,errs,...
%                                                                              SRCavailnum(ipf,:),SRCfixnum2(ipf),SRCsucnum2(ipf),SRCfixnum1(ipf),SRCsucnum1(ipf));
%          end
                                                                           
%         FARmus
    end
end


%% 3. calculate the statistics
nstemp=zeros(Pfttlslen,1);
% ARmethodInd=1;%FAR
% [resFAR]=collectstats(na*ones(Pfttlslen,1),nstemp,FARfixnum1,FARsucnum1,FARfixnum2,FARsucnum2,FARavailnum,Nsamp,ep,Pfttls,Psb,ARmethodInd);
% ARmethodInd=2; %SRC
% [resSRC]=collectstats(nstemp,ns2,SRCfixnum1,SRCsucnum1,SRCfixnum2,SRCsucnum2,SRCavailnum,Nsamp,ep,Pfttls,Psb,ARmethodInd);
ARmethodInd=5; %Rebuttal TSRC benchmark
[resbenchmark]=collectstats(optsub,ns2,fixnum1,sucnum1,fixnum2,sucnum2,availnum,Nsamp,ep,Pfttls,Psb,ARmethodInd);
ARmethodInd=6; %Rebuttal TSRC fit
[resfitTSRC]=collectstats(Fitoptsub,ns2,Fitfixnum1,Fitsucnum1,Fitfixnum2,Fitsucnum2,Fitavailnum,Nsamp,ep,Pfttls,Psb,ARmethodInd);
resallARmethods=[resbenchmark;resfitTSRC];

display(ep);
%%End of the file-0222220**********************9
