function usdt=runonce_ffrt_tsrc(num)
tic
S=load('options.mat');
opts=S.opts;
clear S;
filename=opts(num).filename;
Nsamp=opts(num).Nsamp;
Nsamp=1E6;%debug
loadfile=strcat('./simumat/',filename,'.mat');
if ~exist(loadfile,'file')
    display(strcat('File doesnt exist: ',loadfile));
    usdt=0;
    return;
end
S=load(loadfile);
res=S.res; clear S;
lenres=length(res);
lenres=48;%to shorten the running time.
savefileffrt=strcat('./resFFRT/','FFRT_',filename,'.txt');
% savefiletsrc=strcat('TSRC_',filename,'.txt');
% check if the epoch has been simulated
epochstosimu=[];
for i=1:lenres
    tempfileffrt=strcat('./resFFRT/','FFRTModel',num2str(num),'Epoch',num2str(i),'.txt');
    if ~exist(tempfileffrt,'file')
        epochstosimu=[epochstosimu, i];
    end
end

epochsimulens=length(epochstosimu);
display('initial time');
toc
% P0=0;%0.9;
% parfor i=1:lenres
parfor epi=1:epochsimulens

i=epochstosimu(epi);
tempfileffrt=strcat('./resFFRT/','FFRTModel',num2str(num),'Epoch',num2str(i),'.txt');
% tempfiletsrc=strcat('TSRCModel',num2str(num),'Epoch',num2str(i),'.txt');

% display(i);%FixRateOneEp(ep,Qa,Nsamp,Psb,P0)
[resffrt]=ffrt_OneEp(i,res(i).Qa,Nsamp);
% [resffrt, restsrc]=muFixRateOneEp(i,res(i).Qa,res(i).Qb,res(i).Qab,Nsamp,res(i).Ps,P0);
% [reseprt01,reseprt001,resepdt01,resepdt001]=muFixRateOneEp(i,res(i).Qa,res(i).Qb,res(i).Qab,Nsamp,res(i).Ps,0);

fid1=fopen(tempfileffrt,'w');
fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',resffrt');
fclose(fid1);

% fid2=fopen(tempfiletsrc,'w');
% fprintf(fid2,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',restsrc');
% fclose(fid2);

end


display('after parfor');
toc
slicefiles1=cell(16,1);
% slicefiles2=cell(16,1);

lenslice=lenres/16;
parfor i=1:16
    slicefiles1{i,1}=strcat('./resFFRT/','FFRTSliceModel',num2str(num),'Slice',num2str(i),'.txt');
%     slicefiles2{i,1}=strcat('TSRCSliceModel',num2str(num),'Slice',num2str(i),'.txt');
    
    epcs=lenslice*(i-1)+1:lenslice*i;
    
    ffrt_slicewrite(slicefiles1{i,1},'./resFFRT/FFRT','',num,epcs);
%     tsrcslicewrite(slicefiles2{i,1},'TSRC','',num,epcs);

%     write to 16 slices
end
display('SliceFiles written');
toc

% numslice=min(numslice1,numslice2,numslice3,numslice4);

fid1=fopen(savefileffrt,'w');
% fid2=fopen(savefiletsrc,'w');

for i=1:16
    tempfileffrt=slicefiles1{i,1};
%     tempfiletsrc=slicefiles2{i,1};
%     if ~exist(tempfileffrt)|| ~exist(tempfiletsrc)
    if ~exist(tempfileffrt,'file')
        display(strcat('File doesnt exist: ',tempfileffrt));
        continue;
    end
    reseprt01=load(tempfileffrt);
    fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt01');
    delete(tempfileffrt);
    
%     reseprt001=load(tempfiletsrc);
%     fprintf(fid2,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt001');
%     delete(tempfiletsrc);
    
    %     write to one file
end

fclose(fid1);
% fclose(fid2);
%completefile=strcat('comp',savefileffrt);
%movefile(savefileffrt,completefile);
display('all files combined');
toc
usdt=toc;
