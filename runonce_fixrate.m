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
    display(strcat('File doesnt exist: ',loadfile));
    usdt=0;
    return;
end
S=load(loadfile);
res=S.res; clear S;
lenres=length(res);
savefilert01=strcat('RT_',filename,'Pf01.txt');
savefilert001=strcat('RT_',filename,'Pf001.txt');

savefiledt01=strcat('DT_',filename,'Pf01.txt');
savefiledt001=strcat('DT_',filename,'Pf001.txt');
display('initial time');
toc
P0=0.9;
parfor i=1:lenres
tempfilert01=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
tempfilert001=strcat('RTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');

tempfiledt01=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf01.txt');
tempfiledt001=strcat('DTModel',num2str(num),'Epoch',num2str(i),'Pf001.txt');

% display(i);%FixRateOneEp(ep,Qa,Nsamp,Psb,P0)
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
display('after parfor');
toc
slicefiles1=cell(16,1);
slicefiles2=cell(16,1);
slicefiles3=cell(16,1);
slicefiles4=cell(16,1);

lenslice=lenres/16;
parfor i=1:16
    slicefiles1{i,1}=strcat('RTSliceModel',num2str(num),'Slice',num2str(i),'Pf01.txt');
    slicefiles2{i,1}=strcat('RTSliceModel',num2str(num),'Slice',num2str(i),'Pf001.txt');
    slicefiles3{i,1}=strcat('DTSliceModel',num2str(num),'Slice',num2str(i),'Pf01.txt');
    slicefiles4{i,1}=strcat('DTSliceModel',num2str(num),'Slice',num2str(i),'Pf001.txt');
    
    epcs=lenslice*(i-1)+1:lenslice*i;
    
    slicewrite(slicefiles1{i,1},'RT','Pf01',num,epcs);
    slicewrite(slicefiles2{i,1},'RT','Pf001',num,epcs);
    slicewrite(slicefiles3{i,1},'DT','Pf01',num,epcs);
    slicewrite(slicefiles4{i,1},'DT','Pf001',num,epcs);

%     write to 16 slices
end
display('SliceFiles written');
toc

% numslice=min(numslice1,numslice2,numslice3,numslice4);

fid1=fopen(savefilert01,'w');
fid2=fopen(savefilert001,'w');
fid3=fopen(savefiledt01,'w');
fid4=fopen(savefiledt001,'w');

for i=1:16
    tempfilert01=slicefiles1{i,1};
    tempfilert001=slicefiles2{i,1};
    tempfiledt01=slicefiles3{i,1};
    tempfiledt001=slicefiles4{i,1};
    if ~exist(tempfilert01)|| ~exist(tempfilert001)|| ~exist(tempfiledt01)|| ~exist(tempfiledt001)
        display(strcat('File doesnt exist: ',tempfilert01));
        continue;
    end
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
    %     write to one file
end

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
display('all files combined');
toc
usdt=toc;
