function usdt=distcomb(modelnum,ep)
tic
% cluster1 = parcluster('Iridis4_16node_16h');
% cluster1=parcluster('local');%debug
S=load('options.mat');
opts=S.opts;
clear S;

% lenepochs=1;
totalNsamp=1E9;
% totalNsamp=3000;%debug
% Nblocks=3;

filename=opts(modelnum).filename;
loadfile=strcat('./simumat/',filename,'.mat');
if ~exist(loadfile,'file')
    display(strcat('File doesnt exist: ',loadfile));
    usdt=0;
    return;
end
S=load(loadfile);
res=S.res; clear S;
%change the dir name, and the filenamepre
filenamepre=strcat('./resFFRTE9/','FFRTModel',num2str(modelnum));

filenameepc=strcat(filenamepre,'Epoch',num2str(ep),'.txt');
combmat(filenameepc,ep,res(ep).Qa,totalNsamp);%Qa,Qb,Qab,Ps

%         if ~exist(filenameepc,'file')
% %             [epcusdt]=ffrt_falsealarm_OneEp(filenameepc,ep,res(ep).Qa,totalNsamp);
% %debug
%           job = batch(cluster1,@ffrtRatioDist_parfor_OneEp,1,{filenameepc,ep,res(ep).Qa,totalNsamp,Nblocks},'Pool',3);
%         end

toc
usdt=toc;