function [likSum ] = logSum(liks,weights,isLogWeights)
%% logSum: log(sum(exp(liks)*weights)) for posterior predictives
% liks is normalized to max to preserve numerical accuracy

if nargin<2
    weights = log(ones(size(liks))/size(liks,2)); %equally likely
end
if nargin<3 || isempty(isLogWeights)
    isLogWeights=1;
end
if ~isLogWeights
    weights=log(weights);
end

temp=liks+weights;
maxLik = max(temp(:));
normLiks = temp-maxLik; % normalized to preserve numerical accuracy

%NB - bound for numerical robustness
lmax = log(realmax);
lmin= log(realmin);
maxLik(maxLik<lmin)=lmin; 
maxLik(maxLik>lmax)=lmax;

expLiks = exp(normLiks); 
expAdj = exp(maxLik);
likSum1 = sum(expLiks,2)*expAdj; %rows represent same state
likSum = log(likSum1);