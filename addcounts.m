function [counts]=addcounts(counts,counttmp)

len=length(counts);
for ns=1:len
    counts(ns).sucnumILS=counts(ns).sucnumILS+counttmp(ns).sucnumILS;
    counts(ns).fixnumrt=counts(ns).fixnumrt+counttmp(ns).fixnumrt;
    counts(ns).sucnumrt=counts(ns).sucnumrt+counttmp(ns).sucnumrt;
    counts(ns).failnumrt=counts(ns).failnumrt+counttmp(ns).failnumrt;
    counts(ns).falsealarmrt=counts(ns).falsealarmrt+counttmp(ns).falsealarmrt;
end