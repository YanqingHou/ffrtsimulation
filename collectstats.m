function [resavails]=collectstats(ns1,ns2,fixnum1,sucnum1,fixnum2,sucnum2,availnum,Nsamp,ep,Pfttls,Psb,ARmethodInd)


failnum1=fixnum1-sucnum1;
failnum2=fixnum2-sucnum2;


fixrate1=fixnum1/Nsamp;
sucrate1=sucnum1/Nsamp;
failrate1=failnum1/Nsamp;
fixrate2=fixnum2./(Nsamp-fixnum1);
sucrate2=sucnum2./(Nsamp-fixnum1);
failrate2=failnum2./(Nsamp-fixnum1);

fixratettl=(fixnum1+fixnum2)/Nsamp;
sucratettl=(sucnum1+sucnum2)/Nsamp;
failratettl=(failnum1+failnum2)/Nsamp;
availability=availnum/Nsamp;
% resbenchmark should also contain required total failure rate.
Pfttlslen=length(Pfttls);
resavails=[ep*ones(Pfttlslen,1),ARmethodInd*ones(Pfttlslen,1),Psb*ones(Pfttlslen,1), Pfttls, ns1,ns2,fixrate1, sucrate1, failrate1,fixrate2,sucrate2,failrate2, fixratettl, sucratettl, failratettl, availability];
nanidx=isnan(resavails);
resavails(nanidx)=0;