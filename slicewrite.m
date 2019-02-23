function slicenum=slicewrite(filename,pref,postf,num,epcs)
% slicefiles1{i,1},'RT','Pf01',num,epcs
slicenum=0;
fpw=fopen(filename,'w');
for i=epcs
    tempfile=strcat(pref,'Model',num2str(num),'Epoch',num2str(i),postf,'.txt');
    if ~exist(tempfile)
        display(strcat('File doesnt exist: ',tempfile));
        continue;
    end    
    slicenum=slicenum+1;
    reseprt=load(tempfile);
    fprintf(fpw,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt');
    delete(tempfile);
end    

fclose(fpw);