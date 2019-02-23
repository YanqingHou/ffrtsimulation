function Mu2=cpmu2(Pf1tol,Pf2tol,ns)

if ns<=0 || ns>65 || ~isnumeric(ns)
  error('ns is not valid, or should be positive and less than 65!');
end

S=load('FitFuncsMu2PfILS.mat');
fitfuncs=S.fitfuncs;
clear S;

Pf1s=[0.0009,0.009];
% Pfslen=length(Pfs);

[idxpf,~]=find(Pf1s<=Pf1tol,1,'last');
if isempty(idxpf)
   Mu2=0;
   return
end

if Pf2tol>=0.1 % PfILS>=0.2
    Mu2=fitfuncs{ns,idxpf}(0.1);
elseif Pf2tol>0 % Pftol<=PfILS<0.2
    Mu2=fitfuncs{ns,idxpf}(Pf2tol);  
else % PfILS<Pftol
   Mu2=0;
end

if Mu2>1
    Mu2=1;
elseif Mu2<0
    Mu2=0;
else
end
%Mu2=0.9;
%Mu2=0.95;
