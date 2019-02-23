cluster1 = parcluster ;
[pending, queued, running, completed]=findJob(cluster1);
% completed=findJob(cluster1);
for jobid=1:length(completed)
joboutput = getAllOutputArguments(completed(jobid));
% joboutput = fetchOutputs(job);
savefile=joboutput{1,2};
lenres=size(joboutput,1);%length(jobs(jobid).Tasks);
delete(completed(jobid));
fid=fopen(savefile,'w');
for i=1:lenres
%     tempfile=strcat('Model',num2str(num),'Epoch',num2str(i),'.txt');
     resep=joboutput{i,1};
     fprintf(fid,'%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n',resep');
%     delete(tempfile);
end
fclose(fid);
end