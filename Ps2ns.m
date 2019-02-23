function ns=Ps2ns(Psb,Ps0,na,D)
    k1 = 1;
    while Psb < Ps0 && k1 < na  
        k1 = k1 + 1;
        % bootstrapped success rate if the last n-k+1 ambiguities would be fixed
        Psb = prod ( 2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );   
    end
    if Psb>=Ps0
        ns=na-k1+1;
    else
        ns=0;
    end
    
    return