function batchsimu_tsrc_hours(alljob, hours)
% cluster1 = parcluster ;
%cluster1.SubmitArguments='-l walltime=2:00:00, mem=4gb';
% alljob=[141:160];
switch hours
    case 0.1
        cluster1 = parcluster('Iridis4_1node_5min');
    case 4
        cluster1 = parcluster('Iridis4_1node_4h');
    case 8
        cluster1 = parcluster('Iridis4_1node_8h');
    case 12
        cluster1 = parcluster('Iridis4_1node_12h');
    case 16
        cluster1 = parcluster('Iridis4_1node_16h');
    case 24
        cluster1 = parcluster('Iridis4_1node_24h');
    case 48
        cluster1 = parcluster('Iridis4_1node_48h');
    otherwise
        cluster1 = parcluster;
end
for i=1:length(alljob)
    job = batch(cluster1,@runonce_newtsrc,1,{alljob(i)},'Pool',15);
end
