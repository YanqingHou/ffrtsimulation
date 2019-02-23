function usdt=runonce_tsrc(num)
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
loadfileffrt=strcat('./resFFRT/FFRT_',filename,'.txt');
savefiletsrc=strcat('./resTSRC/TSRC_',filename,'.txt');

display('initial time');
toc


if ~exist(loadfileffrt,'file');
    display(strcat('File doesnt exist: ', loadfileffrt));
    usdt=0;
    return;
end

ffrtclnep=1;
PfixData0=load(loadfileffrt);

slicedPfixData=cell(lenres,1);
for i=1:lenres
    eprows=PfixData0(:,ffrtclnep)==i;
    slicedPfixData{i,1}=PfixData0(eprows,:);
end
clear PfixData0;
parfor i=1:lenres

tempfileffrt=strcat('./resTSRC/TSRCModel',num2str(num),'Epoch',num2str(i),'.txt');
[restsrc]=tsrc_OneEp(i,res(i).Qa,Nsamp,res(i).Ps,slicedPfixData{i,1});

fid1=fopen(tempfileffrt,'w');
fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',restsrc');
fclose(fid1);

end

clear PfixData;
display('after parfor');
toc
slicefiles1=cell(16,1);
% slicefiles2=cell(16,1);

lenslice=lenres/16;
parfor i=1:16
    slicefiles1{i,1}=strcat('./resTSRC/TSRCSliceModel',num2str(num),'Slice',num2str(i),'.txt');
%     slicefiles2{i,1}=strcat('TSRCSliceModel',num2str(num),'Slice',num2str(i),'.txt');
    
    epcs=lenslice*(i-1)+1:lenslice*i;
    
    tsrc_slicewrite(slicefiles1{i,1},'./resTSRC/TSRC','',num,epcs);
%     tsrcslicewrite(slicefiles2{i,1},'TSRC','',num,epcs);

%     write to 16 slices
end
display('SliceFiles written');
toc

% numslice=min(numslice1,numslice2,numslice3,numslice4);

fid1=fopen(savefiletsrc,'w');
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
    fprintf(fid1,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt01');
    delete(tempfileffrt);

    
    %     write to one file
end

fclose(fid1);
% fclose(fid2);
display('all files combined');
toc
usdt=toc;
