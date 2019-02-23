function batchsimu_forpar(alljob)
cluster1 = parcluster ;
%cluster1.SubmitArguments='-l walltime=2:00:00, mem=4gb';
% alljob=[141:160];
for i=1:length(alljob)
    job = batch(cluster1,@runonce_fixrate,1,{alljob(i)},'Pool',15);
end
