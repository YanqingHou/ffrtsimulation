function svmu=svcpmu(Pf_ILS,Pf,ns)
if Pf_ILS<=Pf
    svmu=1;
else
    svmu=ratioinv(Pf,Pf_ILS,ns);
end