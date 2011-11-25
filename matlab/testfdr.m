function k=testfdr(pvals,q)
%Benjamini-Hochberg method for controlling false discovery rate
%SYNTAX
%k=testfdr(pvals,q) 
%INPUTS
%pvals is a column vector of p-values
%q is the level at which we want to control FDR
%OUTPUTS
%k is a column vector with zero for null hypotheses we shouldn't reject and
%1 for null hypotheses we should reject, for expected false discovery rate

[p,ind]=sort(pvals);
m=length(pvals);
i=(1:m)';
imq=i/m*q;
ki=find(p<=imq,1,'last');%LARGEST i satisfying condition
k=zeros(m,1);
k(1:ki)=1;
k(ind)=k;%original order of p-values

