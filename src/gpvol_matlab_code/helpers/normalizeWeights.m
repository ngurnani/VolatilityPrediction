function [normWgts Neff] = normalizeWeights(logDensity,dim)

if nargin<2
    dim=1; %1 by row; 2 by col
end

dataSize=size(logDensity);
repSize = ones(1,length(dataSize));
repSize(dim)=dataSize(dim);

maxLogDensity = max(logDensity,[],dim);
%adjust for numerical stability
wgts = exp(logDensity - repmat(maxLogDensity,repSize )); 

%nan or imaginary wgts is bad support or extreme low probability
wgts(imag(wgts)~=0)=0; 
normWgts = wgts./repmat(nansum(wgts,dim),repSize); 
normWgts(isnan(normWgts))=0;

%track effective sample size    
Neff  = 1./sum(normWgts.^2,dim); 