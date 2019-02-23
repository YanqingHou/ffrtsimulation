function usdt=runonce_bench(num)
%% Initialization
tic
S=load('options.mat');
opts=S.opts;
clear S;
filename=opts(num).filename;
%Nsamp=opts(num).Nsamp;
Nsamp=1E6;%debug
loadfile=strcat('./simumat/',filename,'.mat');
if ~exist(loadfile,'file')
    display(strcat('File doesnt exist: ',loadfile));
    usdt=0;
    return;
end
S=load(loadfile);
res=S.res; clear S;
%lenres=length(res);
lenres=48;%to shorten the running time.
loadfileffrt=strcat('./resFFRT/FFRT_',filename,'.txt');
loadfiletsrc=strcat('./resnewTSRC/TSRC_',filename,'.txt');

savefilebench=strcat('./resBench/Bench_',filename,'.txt');

display('initial time');
toc

%% load existing files
if ~exist(loadfileffrt,'file')
    display(strcat('File doesnt exist: ', loadfileffrt));
    usdt=0;
    return;
end

if ~exist(loadfiletsrc,'file')
    display(strcat('File doesnt exist: ', loadfiletsrc));
    usdt=0;
    return;
end
ffrtclnep=1;

PfixData0=load(loadfileffrt);
slicedPfixData=cell(lenres,1);

TSRCData0=load(loadfiletsrc);
slicedTSRCData=cell(lenres,1);
for i=1:lenres
    eprows=PfixData0(:,ffrtclnep)==i;
    slicedPfixData{i,1}=PfixData0(eprows,:);
    eprows=TSRCData0(:,ffrtclnep)==i;
    slicedTSRCData{i,1}=TSRCData0(eprows,:);
end
clear PfixData0; clear TSRCData0; clear eprows;
%% Parallel running
% parfor i=1:lenres
parfor i=1:lenres

tempfilebench=strcat('./resBench/BenchModel',num2str(num),'Epoch',num2str(i),'.txt');
[resbench]=benchmark_OneEp(i,res(i).Qa,res(i).Qb,res(i).Qab,Nsamp,res(i).Ps,slicedPfixData{i,1},slicedTSRCData{i,1});

fid1=fopen(tempfilebench,'w');
fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',resbench');
fclose(fid1);

end

clear slicedPfixData; clear slicedTSRCData;
display('after parfor');
toc
%% Parallel saving result files
slicefiles1=cell(16,1);
% slicefiles2=cell(16,1);

lenslice=lenres/16;
parfor i=1:16
    slicefiles1{i,1}=strcat('./resBench/BenchSliceModel',num2str(num),'Slice',num2str(i),'.txt');
%     slicefiles2{i,1}=strcat('TSRCSliceModel',num2str(num),'Slice',num2str(i),'.txt');
    
    epcs=lenslice*(i-1)+1:lenslice*i;
    
    bench_slicewrite(slicefiles1{i,1},'./resBench/Bench','',num,epcs);
%     tsrcslicewrite(slicefiles2{i,1},'TSRC','',num,epcs);

%     write to 16 slices
end
display('SliceFiles written');
toc
%% Combine all result files into one
% numslice=min(numslice1,numslice2,numslice3,numslice4);

fid1=fopen(savefilebench,'w');
% fid2=fopen(savefiletsrc,'w');

for i=1:16
    tempfilebench=slicefiles1{i,1};
%     tempfiletsrc=slicefiles2{i,1};
%     if ~exist(tempfileffrt)|| ~exist(tempfiletsrc)
    if ~exist(tempfilebench,'file')
        display(strcat('File doesnt exist: ',tempfilebench));
        continue;
    end
    reseprt01=load(tempfilebench);
    fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt01');
    delete(tempfilebench);

    
    %     write to one file
end

fclose(fid1);
% fclose(fid2);
display('all files combined');
toc
usdt=toc;
