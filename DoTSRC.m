function [availnum,fixnum2,sucnum2,fixnum1,sucnum1]=DoTSRC(zhat,ztrue,bhat,L,D,na,ns1,ns2,mu1,mu2,...
                                                                             Kcoef1,Kcoef2,beflt,errs,...
                                                                             availnum,fixnum2,sucnum2,fixnum1,sucnum1)
ncands=2;
% for ipf=1:Pfttlslen
        if ns1>0
            k1=na-ns1+1;
            [zpar1,sqnorm1] = ssearch(zhat(k1:end),L(k1:end,k1:end),D(k1:end),ncands);
            sucflag1=sum(~(ztrue(k1:end)==zpar1(:,1)))==0;
            if sqnorm1(1)/sqnorm1(2)>mu1
                if ns2>0
                    k2=na-ns2+1;
                    [zpar2,sqnorm2] = ssearch(zhat(k2:end),L(k2:end,k2:end),D(k2:end),ncands);
                    sucflag2=sum(~(ztrue(k2:end)==zpar2(:,1)))==0;
                    %                     ratio test for the second step.
                    if sqnorm2(1)/sqnorm2(2)>mu2
                        %                             use the float solution, update the
                        %                             availability
%                         availnum(ipf,:)=availnum(ipf,:)+double(beflt<=errs);
                        availnum=availnum+double(beflt<=errs); 
                        
                    else
                        fixnum2=fixnum2+1;
                        sucnum2=sucnum2+sucflag2;
%                         failnum2=fixnum2-sucnum2;
                        bchk2=bhat-Kcoef2*(zhat(k2:end)-zpar2(:,1));
                        availnum=availnum+double(norm(bchk2(1:3))<=errs);
                        %                             update the baseline solution with ns2(ipf)
                    end
                else % ns2==0, use float solution.
                    availnum=availnum+double(beflt<=errs);
                end
            else % subset in the first step accepted
                fixnum1=fixnum1+1;
                sucnum1=sucnum1+sucflag1;
%                 failnum1=fixnum1-sucnum1;
                bchk1=bhat-Kcoef1*(zhat(k1:end)-zpar1(:,1));
                availnum=availnum+double(norm(bchk1(1:3))<=errs);
                %                       use the optimal subset
                %                       collect the fix num, success num and fail num
                %                       calculate the availability
            end
            
        else % optsub(ipf)==0
            %                   use the float solution
            availnum=availnum+double(beflt<=errs);
        end