function batchsimu_condpf_ffrt_hours(alljob, hours)
% cluster1 = parcluster ;
%cluster1.SubmitArguments='-l walltime=2:00:00, mem=4gb';
% alljob=[141:160];
switch hours
    case 0.1
        cluster1 = parcluster('Iridis4_1node_5min');
    case 2
	    cluster1 = parcluster('Iridis4_1node_2h');
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

GSF=[1:20,181:200,361:380]; GDF=[21:40,201:220,381:400]; GTF=[41:60,221:240,401:420];
BSF=[61:80,241:260,421:440];BDF=[81:100,261:280,441:460];BTF=[101:120,281:300,461:480];
GBSF=[121:140,301:320,481:500];GBDF=[141:160,321:340,501:520];GBTF=[161:180,341:360,521:540];
strongm=[1:5:16]; mediumm=[2,3,7,8,12,13,17,18]; weakm=[4,5,9,10,14,15,19,20];

if alljob==0
    alljob=[GBDF,GBTF];%[GBDF(mediumm),GBDF(weakm)];
end
for i=1:length(alljob)
    job = batch(cluster1,@runonce_ffrt_condpf,1,{alljob(i)},'Pool',15);
end
