function [uniques] = numUnique(data,dim)
if nargin<2
    dim=2; %across columns be default
end
if dim==1
    data=data'; %execute for columns
end
[t d] = size(data);
uniques=zeros(t,1);
for i=1:t
    uniques(i)=numel(unique(data(i,:)));
end