function ns2=findns2(Psb,P0,na,D)
    k1 = 1;
    while Psb < P0 && k1 < na  
        k1 = k1 + 1;
        % bootstrapped success rate if the last n-k+1 ambiguities would be fixed
        Psb = prod ( 2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );   
    end
    if Psb>=P0
        ns2=na-k1+1;
    else
        ns2=0;
    end
    
    return