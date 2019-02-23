function [usdt]=ffrt_parfor_OneEp(filenameepc,ep,Qa,totalNsamp)%Qa,Qb,Qab,Ps
% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; 
disp('entering new ffrt_parfor_OneEp...');
Nblocks=500;%debug
Nsamp=totalNsamp/Nblocks;
na      = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
% ncands=2;
[Qzhat,~,L,D,ztrue,~] = decorrel(Qa,atrue);
% Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');
counts = repmat(struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[],...
    'themus',[],'thePfix',[],'thePfs',[],'thePscons',[],'theFalseAlarm',[]), na, 1 );

% counttmp=struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[],...
%     'themus',[],'thePfix',[],'thePfs',[],'thePscons',[],'theFalseAlarm',[]);
% counttmp=counts;
% 1. FFRT simulation
% 1.1 Initialization for FFRT
mus=0.00001:0.00001:1;
%mus=0.001:0.001:1;
mus=mus';
muslen=length(mus);

%Pfs=[0.001,0.01]';
Pfs=[1e-7, 1e-5]';%[0.001,0.01]';
% Pfs=[0.0005:0.0001:0.0009,0.001:0.001:0.01]';

Pfslen=length(Pfs);

ffrtclnep=1;ffrtclnna=2;ffrtclnns=3;ffrtclnPsb=4;ffrtclnPsILS=5;ffrtclnPf_req=6;
ffrtclnmu=7;ffrtclnPfix=8;ffrtclnPscon=9; ffrtclnPftrue=10; ffrtclnPfalseAlarm=11;

resffrt=zeros(na*Pfslen,11);
% resffrt=zeros(Pfslen,11);


% if na<8
%     return;
% end
% 1.2 initialization of the struct array
for ns=1:na
    k=na-ns+1;
    counts(ns).Psb=prod(2 * normcdf(1./(2*sqrt(D(k:end)))) -1 );%scalar
    counts(ns).sucnumILS=0;%scalar
   
    counts(ns).fixnumrt=zeros(muslen,1);%vector (muslen x 1)
    counts(ns).sucnumrt=counts(ns).fixnumrt;%vector (muslen x 1)
    counts(ns).failnumrt=counts(ns).fixnumrt;%vector (muslen x 1)
    counts(ns).falsealarmrt=counts(ns).fixnumrt;%vector (muslen x 1)
end
% cellcounts=cell(na,Nblocks);
% ÿ����Ԫ����һ��job,ÿ����Ԫ������500��block������Ҫ500��processors��
% ÿ��block����2*10^6��������Ҫ16Сʱ����ÿ������ڵ�ֻ��16��core��
% ����ÿ��job��Ҫ500/16=31.25=32���ڵ㡣32���ڵ���iridis4�ܹ����뵽�������ޡ�
% ��ʱ������walltime���60Сʱ��walltimeԽС��Խ�����Ŷӣ�����ҲԽ������ϵͳ������
% ÿ��block����2*10^6�������洢������Ҫ���ڴ�Ϊ256MB������ÿ��processor��Ҫ����1G�ڴ档

% ÿ����Ԫ����һ��job,ÿ����Ԫ������250��block������Ҫ250��processors��
% ÿ��block����4*10^6��������Ҫ32Сʱ����ÿ������ڵ�ֻ��16��core��
% ����ÿ��job��Ҫ250/16=15.625=16���ڵ㡣32���ڵ���iridis4�ܹ����뵽�������ޡ�
% ��ʱ������walltime���60Сʱ��walltimeԽС��Խ�����Ŷӣ�����ҲԽ������ϵͳ������
% ÿ��block����4*10^6�������洢������Ҫ���ڴ�Ϊ512MB������ÿ��processor��Ҫ����2G�ڴ档

% -l walltime=16:00:00 -l nodes=16:ppn=16 mem=32GB
filenameblocks=cell(Nblocks,1);
blockinfos=zeros(Nblocks,1);
t_before_parfor=toc;
display(t_before_parfor);
parfor iblock=1:Nblocks  %Nblock=250
    filenameblocks{iblock}=strcat(filenameepc,'block',num2str(iblock),'.mat');
    blockinfos(iblock)=ffrt_OneBlock(filenameblocks{iblock},Qa,Nsamp);%write the simulation result within the function
%     counts=addcounts(counts,counttmp);
end

t_after_parfor=toc;
display(t_after_parfor);
for iblock=1:Nblocks
    if blockinfos(iblock)==0
        disp(strcat(filenameblocks{iblock},'doesnot exist!'));
        continue;
    else
    load(filenameblocks{iblock});
    counts=addcounts(counts,blockcounts);
    delete(filenameblocks{iblock});
    end
end

t_after_comb=toc;
display(t_after_comb);
% parameter for FFRT
for ns=1:na
    counts(ns).sucnumrt=counts(ns).sucnumrt./counts(ns).fixnumrt; % Pscon
    counts(ns).failnumrt=counts(ns).failnumrt/totalNsamp; % Pf
    counts(ns).fixnumrt=counts(ns).fixnumrt/totalNsamp; % Pfix
    counts(ns).falsealarmrt=counts(ns).falsealarmrt/totalNsamp; % Pfa

    % if the fix rate is 0, set failure rate to 1 and conditaional success
    % rate to 0.
    onesidx=counts(ns).fixnumrt==0; 
    counts(ns).failnumrt(onesidx)=1;
    counts(ns).sucnumrt(onesidx)=0;
    for pfi=1:Pfslen
        [rowtmp, ~, mupf]=find(counts(ns).failnumrt<=Pfs(pfi),1,'last'); % find the critical value with a given failure rate
        if isempty(mupf)
%            mupf=0; 
           counts(ns).themus(pfi)=0;%mupf; mu=0 means reject any candidate
           counts(ns).thePfix(pfi)=0;%counts(ns).fixnumrt(row);
           counts(ns).thePscons(pfi)=0;%counts(ns).sucnumrt(row);
           counts(ns).thePfs(pfi)=0;%counts(ns).failnumrt(row);
           counts(ns).theFalseAlarm(pfi)=0;%counts(ns).failnumrt(row);
        else
            counts(ns).themus(pfi)=mus(rowtmp);
            counts(ns).thePfix(pfi)=counts(ns).fixnumrt(rowtmp);
            counts(ns).thePscons(pfi)=counts(ns).sucnumrt(rowtmp);
            counts(ns).thePfs(pfi)=counts(ns).failnumrt(rowtmp);
            counts(ns).theFalseAlarm(pfi)=counts(ns).falsealarmrt(rowtmp);
        end
    end

%     resffrt(ns,:)=[ep,counts(ns).Pftruert(pfi),];

%         indffrt=1:Pfslen;
        indffrt=(ns-1)*Pfslen+1:ns*Pfslen;

        resffrt(indffrt,ffrtclnep)=ep; 
        resffrt(indffrt,ffrtclnna)=na; 

        resffrt(indffrt,ffrtclnns)=ns;
        resffrt(indffrt,ffrtclnPsb)=counts(ns).Psb;
        resffrt(indffrt,ffrtclnPsILS)=counts(ns).sucnumILS/totalNsamp;
        resffrt(indffrt,ffrtclnPf_req)=Pfs;
        resffrt(indffrt,ffrtclnmu)=counts(ns).themus;
        resffrt(indffrt,ffrtclnPfix)=counts(ns).thePfix;
        resffrt(indffrt,ffrtclnPscon)=counts(ns).thePscons;
        
        resffrt(indffrt,ffrtclnPftrue)=counts(ns).thePfs;
        resffrt(indffrt,ffrtclnPfalseAlarm)=counts(ns).theFalseAlarm;
end

% write res for each epoch
% filenameepc=strcat(filenamepre,'_epc',num2str(ep));
printformatstr=makeSformat(size(resffrt,2));
% printformatstrs{epi}=printformatstr;

fid=fopen(filenameepc,'w');
fprintf(fid,printformatstr,resffrt');
fclose(fid);

t_after_wfile=toc;
display(t_after_wfile);
disp('leaving new ffrt_parfor_OneEp...');
usdt=toc;
