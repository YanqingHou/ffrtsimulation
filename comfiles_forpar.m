function comfiles_forpar(alljob)
cluster1 = parcluster ;
% cluster1.SubmitArguments='-l walltime=00:30:00, mem=4gb';
% alljob=[141:160];
for i=1:length(alljob)
    job = batch(cluster1,@combitionres,1,{alljob(i)},'Pool',2);
end
