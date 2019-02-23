function usdt=runonce_ffrt_condpf(num)
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
savefileffrt=strcat('./resFFRTcondpf/','condpf_',filename,'.txt');
% savefiletsrc=strcat('TSRC_',filename,'.txt');
% check if the epoch has been simulated
printformatstrs=cell(lenres,1);
epochstosimu=[];
for i=1:lenres
    tempfileffrt=strcat('./resFFRTcondpf/','FFRTcondpfModel',num2str(num),'Epoch',num2str(i),'.txt');
    if ~exist(tempfileffrt,'file')
        epochstosimu=[epochstosimu, i];
    end
end

epochsimulens=length(epochstosimu);
display('initial time');
toc

parfor epi=1:epochsimulens

i=epochstosimu(epi);
tempfileffrt=strcat('./resFFRTcondpf/','FFRTcondpfModel',num2str(num),'Epoch',num2str(i),'.txt');

[resffrt]=ffrt_condpf_OneEp(i,res(i).Qa,Nsamp);

printformatstr=makeSformat(size(resffrt,2));
printformatstrs{epi}=printformatstr;

fid1=fopen(tempfileffrt,'w');
fprintf(fid1,printformatstr,resffrt');
fclose(fid1);

end

printformatstr=printformatstrs{1};

display('after parfor');
toc
slicefiles1=cell(16,1);

lenslice=lenres/16;
parfor i=1:16
    slicefiles1{i,1}=strcat('./resFFRTcondpf/','FFRTcondpfSliceModel',num2str(num),'Slice',num2str(i),'.txt');   
    epcs=lenslice*(i-1)+1:lenslice*i;
    general_slicewrite(slicefiles1{i,1},'./resFFRTcondpf/FFRTcondpf','',num,epcs,printformatstr);
end
display('SliceFiles written');
toc


fid1=fopen(savefileffrt,'w');

for i=1:16
    tempfileffrt=slicefiles1{i,1};
    if ~exist(tempfileffrt,'file')
        display(strcat('File doesnt exist: ',tempfileffrt));
        continue;
    end
    reseprt01=load(tempfileffrt);
    fprintf(fid1,printformatstr,reseprt01');
    delete(tempfileffrt);    
end

fclose(fid1);
display('all files combined');
toc
usdt=toc;
