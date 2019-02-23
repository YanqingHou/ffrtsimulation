function Pfix=cpPfix1(PfILS,Pftol,ns)

if ns<=0 || ns>66 || ~isnumeric(ns)
  error('ns is not valid, or should be positive and less than 66!');
end

S=load('FitFuncsPfixPfILS.mat');
fitfuncs=S.fitfuncs;
clear S;

Pfs=[0.0005, 0.0006, 0.0007, 0.0008, 0.0009,...
    0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01];
% Pfslen=length(Pfs);

[~,idxpf]=find(Pfs<=Pftol,1,'last');
if isempty(idxpf)
   Pfix=0;
   return
end

if PfILS>=0.2
    Pfix=0;
elseif PfILS>=Pftol
    Pfix=fitfuncs{ns,idxpf}(PfILS);  
else
   Pfix=1;
end

if Pfix>1
    Pfix=1;
end
% Pfix1=FitPfixfuncs{ns,idxPfTSRC(ipf)}(Pfreq1);