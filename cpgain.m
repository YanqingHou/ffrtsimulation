function Rs=cpgain(ns,na,D)

    k1=na-ns+1;
    Ps1= prod (2 * normcdf(1./(2*sqrt(D(k1:end)))) -1 );   
    tempvar(nsi).ns=ns;
    tempvar(nsi).Psb=Ps1;
    Qzpar1 = Qzhat(k1:end,k1:end); Zpar1 = Z(:,k1:end); Qbzpar1=Qab'*Zpar1;
    Kcoef1=Qbzpar1/Qzpar1;  
%     resep(resind1,clnns1)=ns; 
    Qbp1=Qb-Kcoef1*Qbzpar1';
    Rs1=sqrt(det(Qb(1:3,1:3))/det(Qbp1(1:3,1:3)));%large than 1