function cntres=combitionres(num)
S=load('options.mat');
filename=S.opts(num).filename;
clear S;
% savefile=strcat('Res_',filename,'.txt');

loadfile=strcat(filename,'.mat');
S=load(loadfile);
lenres=length(S.res);
clear S;
cntres=0;

savefilert01=strcat('RT_',filename,'Pf01.txt');
savefilert001=strcat('RT_',filename,'Pf001.txt');

savefiledt01=strcat('DT_',filename,'Pf01.txt');
savefiledt001=strcat('DT_',filename,'Pf001.txt');


fid1=fopen(savefilert01,'w');
fid2=fopen(savefilert001,'w');
fid3=fopen(savefiledt01,'w');
fid4=fopen(savefiledt001,'w');
for i=1:lenres
    tempfilert01=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
    tempfilert001=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');
    tempfiledt01=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
    tempfiledt001=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');
    if ~exist(tempfilert01)||~exist(tempfilert001)||~exist(tempfiledt01)||~exist(tempfiledt001)
        continue;
    end
    reseprt01=load(tempfilert01);
    fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f\n',reseprt01');
    delete(tempfilert01);
    cntres=cntres+1;

    reseprt001=load(tempfilert001);
    fprintf(fid2,'%12.5f %12.5f %12.5f %12.5f\n',reseprt001');
    delete(tempfilert001);
    
    resepdt01=load(tempfiledt01);
    fprintf(fid3,'%12.5f %12.5f %12.5f %12.5f\n',resepdt01');
    delete(tempfiledt01);
    
    resepdt001=load(tempfiledt001);
    fprintf(fid4,'%12.5f %12.5f %12.5f %12.5f\n',resepdt001');
    delete(tempfiledt001);
    
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);

% display(strcat(savefile,'is ready.'));
