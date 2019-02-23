function [availnum,fixnum2,sucnum2,fixnum1,sucnum1]=DoSRC(zhat,ztrue,bhat,L,D,na,ns2,...
                                                                             Kcoef2,beflt,errs,...
                                                                             availnum,fixnum2,sucnum2,fixnum1,sucnum1)
ncands=2;
if ns2>0
    k2=na-ns2+1;
   [zpar2,~] = ssearch(zhat(k2:end),L(k2:end,k2:end),D(k2:end),ncands);
   sucflag2=sum(~(ztrue(k2:end)==zpar2(:,1)))==0;
   fixnum2=fixnum2+1;
   sucnum2=sucnum2+sucflag2;
   bchk2=bhat-Kcoef2*(zhat(k2:end)-zpar2(:,1));
   availnum=availnum+double(norm(bchk2(1:3))<=errs);
else
   availnum=availnum+double(beflt<=errs);  
end
    