function [there, notthere]=checkfiles(array)
S=load('options.mat');
opts=S.opts;
clear S;
there=[];
notthere=[];
cntt=0;cntnt=0;
for num=array
filename=opts(num).filename;
savefilert01=strcat('TSRC_',filename,'.txt');
if ~exist(savefilert01,'file')
    notthere=[notthere,num];
    cntnt=cntnt+1;
else
    there=[there, num];
    cntt=cntt+1;
   % error(strcat('File doesnt exist: ',loadfile));
end

end
distxt=[num2str(cntt),' files finished simulation!'];
display(distxt);
distxt=[num2str(cntnt), ' files are not finished!'];
display(distxt);
