function slicenum=ffrt_slicewrite(filename,pref,postf,num,epcs)
% slicewrite for the format in ffrt_OneEp: 10 column data
slicenum=0;
fpw=fopen(filename,'w');
for i=epcs
    tempfile=strcat(pref,'Model',num2str(num),'Epoch',num2str(i),postf,'.txt');
    if ~exist(tempfile,'file')
        display(strcat('File doesnt exist: ',tempfile));
        continue;
    end    
    slicenum=slicenum+1;
    reseprt=load(tempfile);
    fprintf(fpw,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',reseprt');
    delete(tempfile);
end    

fclose(fpw);