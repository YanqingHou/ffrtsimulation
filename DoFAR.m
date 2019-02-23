function [availnum,fixnum2,sucnum2,fixnum1,sucnum1]=DoFAR(zhat,ztrue,bhat,L,D,...
                                                                             mu,Kfar,beflt,errs,...
                                                                             availnum,fixnum2,sucnum2,fixnum1,sucnum1)
ncands=2;
[zfixed,sqnorm] = ssearch(zhat,L,D,ncands);
sucflag=sum(~(ztrue==zfixed(:,1)))==0;

  
if sqnorm(1)/sqnorm(2)<=mu % FFRT passed
   bchk=bhat-Kfar*(zhat-zfixed(:,1));
   availnum=availnum+double(norm(bchk(1:3))<=errs);
   
   fixnum1=fixnum1+1;
   sucnum1=sucnum1+sucflag;
else
   availnum=availnum+double(beflt<=errs); 
end


    