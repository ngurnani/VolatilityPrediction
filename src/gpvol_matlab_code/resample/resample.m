function [parentIndex] = resample(normWgts,option, varargin)
if nargin<2
    option=1; %systematic
end

switch option
    case 1 %systematic
        parentIndex = resampleSystematic(normWgts,varargin{:});
    case 2 %multinomial
        parentIndex = resampleMultinomial(normWgts,varargin{:});
    case 3 %residual 
        parentIndex = resampleResidual(normWgts,varargin{:});
    case 4 %stratified 
        parentIndex = resampleStratified(normWgts,varargin{:});
end