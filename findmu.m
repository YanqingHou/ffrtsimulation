function [mu,Pftrue]=findmu(ns,Pfreq1,PfixData)
% ffrtclnep=1;ffrtclnna=2;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;
% ffrtclnPfix=8;ffrtclnPscon=9; 

ffrtclnns=3;
ffrtclnmu=7;
ffrtclnPftrue=10;
% eprows=PfixData(:,ffrtclnep)==ep;
% PfixData=PfixData(eprows,:);
nsrows=PfixData(:,ffrtclnns)==ns;
PfixData=PfixData(nsrows,:);
[row,~,~]=find(PfixData(:,ffrtclnPftrue)<=Pfreq1,1,'last');
if isempty(row)
    mu=0;
    Pftrue=0;
else
    mu=PfixData(row,ffrtclnmu);
    Pftrue=PfixData(row,ffrtclnPftrue);
end
