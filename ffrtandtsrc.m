function [resffrt,restsrc]=ffrtandtsrc(ep,Qa,Qb,Qab,Nsamp,Psb,P0)
% 先不做FFDT，只做FFRT和TSRC
% 结构
% 1.初始化
% 2.生成浮点解 & LAMBDA 并保存结果
% 3.分拣 FFRT的结果
% 4.分拣 TSRC的结果

% 1. Initialization
% input: Qa, Qb, Qab, Ps
%%%%%%%%%%%%%%%%Initialize the float ambiguities%%%%%%%%%%%%%%%%%%%%%%%%%%%
na      = size(Qa,1);
nb      = size(Qb,1);
ncands=2;

Qx=[Qa,Qab;
    Qab', Qb];
Qx      = (tril(Qx,0)+tril(Qx,-1)');

atrue=randi([-200 200],na,1);% generate random integers
btrue=zeros(nb,1);
xtrue=[atrue;btrue];

[Qzhat,Z,L,D,ztrue,~] = decorrel(Qa,atrue);
Qzhat      = (tril(Qzhat,0)+tril(Qzhat,-1)');

n0=na;% it is used for TSRC
if P0>0
    n0=Ps2ns(Psb,P0,na,D);
end
% 2. 生成浮点解 & LAMBDA 并保存结果
for i=1:Nsamp
%     2.1 生成浮点解
        xhat=mvnrnd(xtrue,Qx,1)';
        ahat=xhat(1:na); bhat=xhat(na+1:na+nb);
        zhat=Z'*ahat;
%     2.2 对每一个子集做LAMBDA
    for ns=1:na
        k=na-ns+1; %zp=zhat(na-ns+1:na);
        [zpar,sqnorm] = ssearch(zhat(k:end),L(k:end,k:end),D(k:end),ncands);
        zpar=zpar(:,1);
        bchk=bhat-counts(ns).Kcoef*(zhat(k:end)-zpar);
        
        sucflag=sum(~(ztrue(k:end)==zpar))==0;                       % 
        counts(ns).sucnumILS = counts(ns).sucnumILS+sucflag;         % Correctly fixed
%         sucflags(ns)=sucflag;
%         counts(ns).fltavail = counts(ns).fltavail+double(bert(1)<=counts(ns).errs');
        counts(ns).ratio=sqnorm(1)/sqnorm(2);
%     2.3 对每个子集做ratio test， 使用100个不同的critical value

        rtpass=sqnorm(1)/sqnorm(2)<counts(ns).rtmus;                % 100 different mu
        counts(ns).fixnumrt=counts(ns).fixnumrt+rtpass;             % accepted by FFRT
        counts(ns).sucnumrt=counts(ns).sucnumrt + double(rtpass&sucflag);   % success rate after FFRT: (Accepted by FFRT & correct) / N
        counts(ns).failnumrt=counts(ns).fixnumrt-counts(ns).sucnumrt;% failure rate after FFRT (Fix-Fix & correct)
%         bert(rtpass,ns)=norm(bchk(1:3));                                % if FFRT passed, baseline estimation is updated, and the rms is updated.
        resmat=comparevect(bert(:,ns),counts(ns).errs);                    % compare the rms with required errors for each mu in FFRT, and this result in a matrix.
        counts(ns).availsrt=counts(ns).availsrt+resmat;                 % In the end, it needs to be devided by N---- the total sample number.
     end
end