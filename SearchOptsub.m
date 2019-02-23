function [optsub,DoStep2flag,paras]=SearchOptsub(Pfttl,Pfreq1,Qzhat,Z,D,Qab,Qb,ns2,g2,PfixData,TSRCData,EstPfixMethod)
%EstPfixMethod: 1--Fitting function; 0--using direct simulated values.

%     find the optimal subset that maximizes Pfix*(g1-g2)
paras=struct('mu1',[],'mu2',[],'K1',[],'K2',[]);
maxdEg=0;
optsub=ns2;
na=length(D);
% initialize K1 as K2. In the optimal subset search process, give it the correct value.
if ns2>0
    k1=na-ns2+1;
    Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
    Kcoef1=Qbzpar1/Qzpar1;
    K2=Kcoef1;
    paras.K2=K2;
    DoStep2flag=1;
    Psb2=prod(2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
    
else
    paras.K2=0;
    DoStep2flag=0;
    Psb2=0;
end
OptK1=paras.K2;
OptPsb= Psb2;
OptPfix1=0;
for ns=ns2+1:na
    k1=na-ns+1;
    
    Psb = prod(2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );
    Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
    Kcoef1=Qbzpar1/Qzpar1;
    Qbp1=Qb-Kcoef1*Qbzpar1';
    g1=sqrt(det(Qb(1:3,1:3))/det(Qbp1(1:3,1:3)))^(1/3);%large than 1
    if EstPfixMethod
        Pfix1=cpPfix1(1-Psb,Pfreq1,ns);
    else
        Pfix1=findPfix(ns,Pfreq1,PfixData);%%%%%%%%%%%%%%
    end
    
    dEg=Pfix1*(g1-g2);
    
    
    if dEg>maxdEg
        maxdEg=dEg;
        optsub=ns;
        OptK1=Kcoef1;
        OptPfix1=Pfix1;
        OptPsb=Psb;
    end
end



paras.K1=OptK1;

if ns2>0          % A1
    %     DoStep2flag=1;   % TSRC
    Pf2req=(Pfttl-Pfreq1)/(1-OptPfix1);
   
    if OptPfix1==1
       paras.mu2=1;
    else 
        if EstPfixMethod
           paras.mu2=cpmu2PfB(Pfreq1,Pf2req,1-Psb2,ns2);
        else
           paras.mu2=findmu2(Pfttl,optsub,ns2,Pf2req,TSRCData);%
        end
    end
    if EstPfixMethod
        paras.mu1=cpmu1(1-OptPsb,Pfreq1,optsub);
    else
        [paras.mu1, ~]=findmu(optsub,Pfreq1,PfixData);
    end
    %     k1=na-ns2+1;
    %     Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
    %     Kcoef1=Qbzpar1/Qzpar1;
    %     K2=Kcoef1;
    %     paras.K2=K2;
    
else %A2
    %     DoStep2flag=0; % ns1 then float
    %     paras.K2=0;
     paras.mu2=0;
     if optsub~=0
        if EstPfixMethod
            paras.mu1=cpmu1(1-OptPsb,Pfreq1,optsub);
        else
            [paras.mu1, ~]=findmu(optsub,Pfreq1,PfixData);
        end
    else
        paras.mu1=0;
    end
    
end
