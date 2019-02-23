function [info]=ffrt_OneBlock(filenameblock,Qa,Nsamp)%Qa,Qb,Qab,Ps
% 0. General Initialization
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
disp('entering new ffrt_OneBlock...');
info=0;
na      = size(Qa,1);
atrue=randi([-200 200],na,1);% generate random integers
ncands=2;
[Qzhat,~,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');

blockcounts = repmat(struct('Psb',[],'ratio',[],'sucnumILS',[],'fixnumrt',[],'sucnumrt',[],'failnumrt',[],'falsealarmrt',[],...
    'themus',[],'thePfix',[],'thePfs',[],'thePscons',[],'theFalseAlarm',[]), na, 1 );
% 1. FFRT simulation
% 1.1 Initialization for FFRT
mus=0.00001:0.00001:1;
%mus=0.001:0.001:1;
mus=mus';
muslen=length(mus);


% 1.2 initialization of the struct array
for ns=1:na
    k=na-ns+1;
    blockcounts(ns).Psb=prod(2 * normcdf(1./(2*sqrt(D(k:end)))) -1 );
    blockcounts(ns).sucnumILS=0;
   
    blockcounts(ns).fixnumrt=zeros(muslen,1);
    blockcounts(ns).sucnumrt=blockcounts(ns).fixnumrt;
    blockcounts(ns).failnumrt=blockcounts(ns).fixnumrt;
    blockcounts(ns).falsealarmrt=blockcounts(ns).fixnumrt;
end

fileprogressname=strcat(filenameblock,'progress.txt');
fid=fopen(fileprogressname,'w');
% 1.3 Simulate all the samples and do LAMBDA, ratio test.
    zhats=mvnrnd(ztrue,Qzhat,Nsamp)';%10^7 samples, used up 1.28GB memory
for i=1:Nsamp
    if ~mod(i,1E4)
    fprintf(fid,'\r%d',i);
    end
    zhat=zhats(:,i);
%     zhat=mvnrnd(ztrue,Qzhat,1)';
    for ns=1:na
        k=na-ns+1; 
        [zpar,sqnorm] = ssearch(zhat(k:end),L(k:end,k:end),D(k:end),ncands);
        zpar=zpar(:,1);
        
        sucflag=sum(~(ztrue(k:end)==zpar))==0;
        blockcounts(ns).sucnumILS = blockcounts(ns).sucnumILS+sucflag;         % Correctly fixed

        blockcounts(ns).ratio=sqnorm(1)/sqnorm(2);
        
        rtpass=blockcounts(ns).ratio<mus;                % 100 different mu
        blockcounts(ns).fixnumrt=blockcounts(ns).fixnumrt+rtpass;             % accepted by FFRT
        blockcounts(ns).sucnumrt=blockcounts(ns).sucnumrt + double(rtpass&sucflag);   % success rate after FFRT: (Accepted by FFRT & correct) / N
        blockcounts(ns).failnumrt=blockcounts(ns).fixnumrt-blockcounts(ns).sucnumrt;% failure rate after FFRT (Fix-Fix & correct)
        blockcounts(ns).falsealarmrt=blockcounts(ns).falsealarmrt + double((~rtpass)&sucflag);%-counts(ns).sucnumrt;% false alarm rate after FFRT (Reject & correct)

    end
 
end
save(filenameblock,'blockcounts');
fclose(fid);
delete(fileprogressname);
t_lv_block=toc;
display(t_lv_block);
disp('leaving new ffrt_OneBlock...');
info=1;
