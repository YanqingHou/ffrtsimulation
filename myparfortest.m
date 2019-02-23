function [consumt,commu]=myparfortest(N)
tic
ticBytes(gcp);
A = 500;
a = zeros(N);
parfor i=1:N
    a(i) = max(abs(eig(rand(A))));
end
consumt=toc
commu=tocBytes(gcp)