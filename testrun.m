
function testrun(N)
% cluster1 = parcluster ;
%cluster1.SubmitArguments='-l walltime=2:00:00, mem=4gb';
% alljob=[141:160];
% switch hours
%     case 0.1
%         cluster1 = parcluster('Iridis4_1node_5min');
%     case 2
% 	    cluster1 = parcluster('Iridis4_1node_2h');
%     case 4
%         cluster1 = parcluster('Iridis4_1node_4h');
%     case 8
%         cluster1 = parcluster('Iridis4_1node_8h');
%     case 12
%         cluster1 = parcluster('Iridis4_1node_12h');
%     case 16
%         cluster1 = parcluster('Iridis4_1node_16h');
%     case 24
%         cluster1 = parcluster('Iridis4_1node_24h');
%     case 48
%         cluster1 = parcluster('Iridis4_1node_48h');
%     otherwise
%         cluster1 = parcluster;
% end
cluster1=parcluster('Iridis4_1node_5min');
N = 200;
% cluster1.SubmitArguments='-l walltime=2:00:00 -l nodes=1:ppn=16 mem=4gb';
cluster1.SubmitArguments = '-l nodes=2:ppn=16,walltime=02:00:00';
% for i=1:length(alljob)
    job = batch(cluster1,@myparfortest,1,{N},'Pool',31);
% end



% [consumt,commu]=myparfortest(N);
