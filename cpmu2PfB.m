function Mu2=cpmu2PfB(Pf1tol,Pf2tol,Pf2B,ns)

if ns<=0 || ns>65 || ~isnumeric(ns)
  error('ns is not valid, or should be positive and less than 65!');
end

S=load('FitFuncsMu2PfBpoly2.mat');

fitfuncs=S.fitfuncs;
clear S;

Pf1s=[0.0009,0.009];
Pf2s=[0.0001:0.0001:0.001,0.002:0.001:0.01,0.02:0.01:0.1]';
% Pf2slen=length(Pf2s);
% Pfslen=length(Pfs);
if Pf2tol>=1
   Mu2=1;
   return; 
end
if Pf2B==0
   Mu2=1;
   return;
end

[idxpftol1,val1]=find(Pf1s<=Pf1tol,1,'last');
if isempty(idxpftol1)
   Mu2=0;
   return
end

[idxpftol2,val2]=find(Pf2s<=Pf2tol,1,'last');
if isempty(idxpftol2)
   Mu2=0;
   return;
else
    fitfun=fitfuncs{idxpftol2,ns,idxpftol1};
    if isnumeric(fitfun)
        Mu2=fitfun;
    else
        Mu2=fitfun(Pf2B);
    end
end

% if Pf2tol>=0.1 % PfILS>=0.2
%     Mu2=fitfuncs{ns,idxpftol1}(0.1);
% elseif Pf2tol>0 % Pftol<=PfILS<0.2
%     Mu2=fitfuncs{ns,idxpftol1}(Pf2tol);  
% else % PfILS<Pftol
%    Mu2=0;
% end

if Mu2>1
    Mu2=1;
elseif Mu2<0
    Mu2=0;
else
end
% Mu2=0.9;
