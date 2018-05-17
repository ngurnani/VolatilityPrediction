function [m ] = rapfPointEstimate(prevParams,prevDensity,lambda)
[N numParams]=size(prevParams);

temp=repmat(exp(prevDensity),1,numParams).*prevParams;
paramMean = sum(temp)/N;

%kernelized smoothing for shared parameters
m = lambda*prevParams+(1-lambda)*repmat(paramMean,N,1); 