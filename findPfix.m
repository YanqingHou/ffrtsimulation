function Pfix1=findPfix(ns,Pfreq1,PfixData)
% ffrtclnep=1;ffrtclnna=2;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;ffrtclnPscon=9;
ffrtclnPfix=8; 
ffrtclnns=3;
% ffrtclnmu=7;
ffrtclnPftrue=10;
% eprows=PfixData(:,ffrtclnep)==ep;
% PfixData=PfixData(eprows,:);
nsrows=PfixData(:,ffrtclnns)==ns;
PfixData=PfixData(nsrows,:);
clear nsrows;
if isempty(PfixData)
   display('In function findPfix ns is not valid!');
   Pfix1=0;
   return;
end
[row,~,~]=find(PfixData(:,ffrtclnPftrue)<=Pfreq1,1,'last');
if isempty(row)
   Pfix1=0; 
else
   Pfix1=PfixData(row,ffrtclnPfix);
end