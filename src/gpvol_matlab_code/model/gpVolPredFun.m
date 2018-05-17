function [vOut] = gpVolPredFun(v,y,gpStruct,hyps,isT)
% gp f on v; v in (-inf, inf)
% v_t = f(h(v_{1:t-1},y_{1:t-1})); some design h of observed and hidden data
% y_t ~ N(0, g(v_{t})); % variance=g(v_{t})
% if v_t=log(sigma_{t}^2) -> g(x) = exp(x)
% if v_t=+-sigma_{t} -> g(x) = x^2

[T N]=size(v);

if nargin<3 || isempty(gpStruct)
    gpStruct.covfunc = {@covMaterniso, 3}; 
    ell = 1/4; sf = 1; 
    gpStruct.hyp.cov = log([ell; sf]);
    
    gpStruct.meanfunc = {@meanSum, {@meanLinear, @meanConst}}; 
    gpStruct.hyp.mean = [0.95; -.1; .1];
    
    gpStruct.likfunc = @likGauss; 
    sn = .1; 
    gpStruct.hyp.lik = log(sn);
    
    hyps = repmat(hyp2Mat(gpStruct.hyp,numel(gpStruct.hyp.mean),numel(gpStruct.hyp.cov)),numel(gpStruct.hyp.lik),[N,1]);
end
meanfunc = gpStruct.meanfunc;
covfunc = gpStruct.covfunc;
likfunc = gpStruct.likfunc;
hyp =gpStruct.hyp; %for formatting purposes

if nargin<5 || isempty(isT)
    isT=0; 
end

vOut=zeros(N,1);
for i=1:N
    curHyp = mat2Hyp(hyps(i,1:end-isT,:),numel(hyp.mean),numel(hyp.cov),numel(hyp.lik));
    
    %construct gp inputs and targets
    [input target]= gpVolInput(v(1:T,i),y(1:T,:)); 
    input=input(1:end-1,:); %gp training input
    target = target(2:end,1); %gp training target
    z = input(end,:); %gp test input
    
    % gp prediction
    [m s2] = gp(curHyp, @infExact, meanfunc, covfunc, likfunc,input,target, z);

    % sample from gp posterior distribution 
    vOut(i)= randn(1)*sqrt(s2)+m; % this allows more particle diversity    
end