function Mu2=cpmu2median(Pf1tol,Pf2tol)


Pf1s=[0.0009,0.009];
% Pfslen=length(Pfs);
coefs=[0.4345, 0.4108, 0.8138;
       0.9434, 0.1063, 0.2523];
[~,idxpf]=find(Pf1s<=Pf1tol,1,'last');
if isempty(idxpf)
   Mu2=0;
   return
end


if Pf2tol<0.0001 
    Pf2tol=0.0001;   
end
 Mu2=coefs(idxpf,1)*Pf2tol^coefs(idxpf,2)+coefs(idxpf,3);

if Mu2>1
    Mu2=1;
elseif Mu2<0
    Mu2=0;
else
    
end
