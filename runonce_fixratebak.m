function usdt=runonce_fixrate(num)
tic
S=load('options.mat');
opts=S.opts;
clear S;
filename=opts(num).filename;
Nsamp=opts(num).Nsamp;
Nsamp=5;%debug
loadfile=strcat(filename,'.mat');
if ~exist(loadfile)
    error(strcat('File doesnt exist: ',loadfile));
end
S=load(loadfile);
res=S.res; clear S;
lenres=length(res);
savefilert01=strcat('RT_',filename,'Pf01.txt');
savefilert001=strcat('RT_',filename,'Pf001.txt');

savefiledt01=strcat('DT_',filename,'Pf01.txt');
savefiledt001=strcat('DT_',filename,'Pf001.txt');

P0=0.9;
parfor i=1:lenres
tempfilert01=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
tempfilert001=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');

tempfiledt01=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
tempfiledt001=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');

display(i);%FixRateOneEp(ep,Qa,Nsamp,Psb,P0)
[reseprt01,reseprt001,resepdt01,resepdt001]=FixRateOneEp(i,res(i).Qa,Nsamp,res(i).Ps,P0);
fid1=fopen(tempfilert01,'w');
fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt01');
fclose(fid1);

fid2=fopen(tempfilert001,'w');
fprintf(fid2,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt001');
fclose(fid2);

fid3=fopen(tempfiledt01,'w');
fprintf(fid3,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',resepdt01');
fclose(fid3);

fid4=fopen(tempfiledt001,'w');
fprintf(fid4,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',resepdt001');
fclose(fid4);
end
% delete(poolobj);
% matlabpool close
fid1=fopen(savefilert01,'w');
fid2=fopen(savefilert001,'w');
fid3=fopen(savefiledt01,'w');
fid4=fopen(savefiledt001,'w');
for i=1:lenres
    tempfilert01=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
    tempfilert001=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');
    tempfiledt01=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
    tempfiledt001=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');
    
    reseprt01=load(tempfilert01);
    fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt01');
    delete(tempfilert01);
    
    reseprt001=load(tempfilert001);
    fprintf(fid2,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt001');
    delete(tempfilert001);
    
    resepdt01=load(tempfiledt01);
    fprintf(fid3,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',resepdt01');
    delete(tempfiledt01);
    
    resepdt001=load(tempfiledt001);
    fprintf(fid4,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',resepdt001');
    delete(tempfiledt001);
    
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
usdt=toc;
