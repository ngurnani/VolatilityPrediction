function[inputs target]=gpVolInput(Vs,data,sigTransform,isSymmetric)
% can add complexity by allowing additional input transformations
if nargin<4 || isempty(isSymmetric)
    isSymmetric=0;
end
if nargin<3 || isempty(sigTransform)
    sigTransform=@(x) x; % identity
end

if isSymmetric
    inputs = [sigTransform(Vs), data.^2];
else
    inputs = [sigTransform(Vs), data];
end
target = sigTransform(Vs);