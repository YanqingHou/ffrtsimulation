function resmat=comparevect(be,errs)
belen=length(be);
errslen=length(errs);
resmat=zeros(belen,errslen);
for i=1:errslen
    resmat(:,i)=be<=errs(i);
end