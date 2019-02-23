function mu2=findmu2(Pfttl,ns1, ns2,Pf2req, TSRCData)
%find mu2 in TSRC result files. make sure the columns are the same as in
%TSRC file

 tsrcclnep=1;tsrcclnna=2;tsrcclnPsb=3;tsrcclnPfttl=4;tsrcclnns=5;tsrcclnPf1=6;tsrcclnPfix1=7;
 tsrcclnns2=8;tsrcclnPssub2=9;tsrcclnPf2req=10;tsrcclnmu2=11;tsrcclnPfix2=12; tsrcclnPf2=13;

ttlrows=TSRCData(:,tsrcclnPfttl)==Pfttl;
TSRCData=TSRCData(ttlrows,:);
if isempty(TSRCData)
   mu2=0;
   display('In function findmu2: Pfttl is not a valid value!');
   return;
end
clear ttlrows;
nsrows=TSRCData(:,tsrcclnns)==ns1;
TSRCData=TSRCData(nsrows,:);
if isempty(TSRCData)
   mu2=0;
   display('In function findmu2: ns1 is not a valid value!');
   return;
end
clear nsrows;
nsrows=TSRCData(:,tsrcclnns2)==ns2;
TSRCData=TSRCData(nsrows,:);
if isempty(TSRCData)
   mu2=0;
   display('In function findmu2: ns2 is not a valid value!');
   return;
end

clear nsrows;

[row,~,~]=find(TSRCData(:,tsrcclnPf2)<=Pf2req,1,'last');
if isempty(row)
    mu2=0;
%     display('cannot find');
else
    mu2=TSRCData(row,tsrcclnmu2);
end