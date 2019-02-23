function [availnum,fixnum2,sucnum2,fixnum1,sucnum1]=DoTSRCstep1(zhat,ztrue,bhat,L,D,na,ns1,ns2,mu1,mu2,...
    Kcoef1,Kcoef2,beflt,errs,...
    availnum,fixnum2,sucnum2,fixnum1,sucnum1)
ncands=2;
if ns1>0
    k1=na-ns1+1;
    [zpar1,sqnorm1] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);
    sucflag1=sum(~(ztrue(k1:end)==zpar1(:,1)))==0;
    if sqnorm1(1)/sqnorm1(2)>mu1
       %                   use the float solution
        availnum=availnum+double(beflt<=errs);
    else % subset in the first step accepted
        fixnum1=fixnum1+1;
        sucnum1=sucnum1+sucflag1;
        %                 failnum1=fixnum1-sucnum1;
        bchk1=bhat-Kcoef1*(zhat(k1:end)-zpar1(:,1));
        availnum=availnum+double(norm(bchk1(1:3))<=errs);
    end
    
else % optsub(ipf)==0
    %                   use the float solution
    availnum=availnum+double(beflt<=errs);
end