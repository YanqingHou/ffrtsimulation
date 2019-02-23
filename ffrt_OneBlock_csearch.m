function [t_lv_block]=ffrt_OneBlock_csearch(filenameblock,Qa,Nsamp)%Qa,Qb,Qab,Ps
% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
disp('entering new ffrt_OneBlock...');
% info=0;
na      = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
ncands=2;
[Qzhat,~,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');

% blockcounts = repmat(struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[],...
%     'themus',[],'thePfix',[],'thePfs',[],'thePscons',[],'theFalseAlarm',[]), na, 1 );
% 1. FFRT simulation
% 1.1 Initialization for FFRT
% mus=0.00001:0.00001:1;
mus=0.001:0.001:1;
mus=mus';
muslen=length(mus);

Psbs=zeros(na,1);
SucnumILSs=Psbs;
fixnumrts=zeros(muslen,na);
sucnumrts=fixnumrts;
% failnumrts=fixnumrts;
falsealarmrts=fixnumrts;
% ratios=nan(na,1);

% 1.2 initialization of the struct array
for ns=1:na
    k=na-ns+1;
    Psbs(ns)=prod(2 * normcdf(1./(2*sqrt(D(k:end)))) -1 );
end

% fileprogressname=strcat(filenameblock,'progress.txt');
% fid=fopen(fileprogressname,'w');
% 1.3 Simulate all the samples and do LAMBDA, ratio test.
%     zhats=mvnrnd(ztrue,Qzhat,Nsamp)';%10^7 samples, used up 1.28GB memory
% load('myzhats.mat');
for i=1:Nsamp
%     if ~mod(i,1E5)
%         str=strcat('In ',filenameblock,' Progress ',num2str(i));
%         disp(str);
%     end
    %     zhat=zhats(:,i);
    zhat=mvnrnd(ztrue,Qzhat,1)';
    for ns=1:na
        k=na-ns+1;
        [zpar,sqnorm] = csearchL(zhat(k:end),L(k:end,k:end),D(k:end),ncands);
        zpar=zpar(:,1);
        
        sucflag=sum(~(ztrue(k:end)==zpar))==0;
        
        SucnumILSs(ns) = SucnumILSs(ns)+sucflag;         % Correctly fixed
        ratio=sqnorm(1)/sqnorm(2);
        rtpass=ratio<mus;
        
        fixnumrts(:,ns)=fixnumrts(:,ns)+rtpass;          % fixed, but correctness unknown
        sucnumrts(:,ns)=sucnumrts(:,ns)+ double(rtpass&sucflag); % fixed and correct
        falsealarmrts(:,ns)=falsealarmrts(:,ns)+ double((~rtpass)&sucflag); % rejected but correct
    end
    
end
failnumrts=fixnumrts-sucnumrts;
blockcounts = repmat(struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[],...
    'themus',[],'thePfix',[],'thePfs',[],'thePscons',[],'theFalseAlarm',[]), na, 1 );

for ns=1:na
    blockcounts(ns).Psb=Psbs(ns);
    blockcounts(ns).sucnumILS=SucnumILSs(ns);
    blockcounts(ns).fixnumrt=fixnumrts(:,ns);
    blockcounts(ns).sucnumrt=sucnumrts(:,ns);
    blockcounts(ns).falsealarmrt=falsealarmrts(:,ns);
    blockcounts(ns).failnumrt=failnumrts(:,ns);
end
save(filenameblock,'blockcounts');
% fclose(fid);
% delete(fileprogressname);
t_lv_block=toc;
display(t_lv_block);
disp('leaving new ffrt_OneBlock...');
% info=1;
