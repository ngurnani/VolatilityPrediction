function [parts paramParts predLik stats] = racpf(data,N,modelStruct,gpStruct,initVs,varargin)
%% input arguments
% data - test data
% N - number of particles to use
% modelStruct - struct of model functions and parameters.
%             - should contain
%                - obsfunc (observation model
%                - transfunc (for manipulating latent input structure)
%                - for GP-Vol parts are natural scale (-inf,inf)
% gpStruct - struct to specify GP functions and hyperparameters. 
%          - should contain
% initVs - initial chain of latent variables 
%% outputs
% parts - particle representation of latent variables
% paramParts - particle represenation of hyper-parameters
% predLik - one-step forward predictive likelihood
% stats - rapcf stats such as:
%           - timeUsed per step
%           - number of effective particles per step. to check degradation of pf
%           - quantiles for the latent variables

[T d] = size(data); %T days, d dimension
parts= zeros(T,N);  %is N*T*d, but d=1 for garch

% how much initial data used to seed the inference.
burnIn=modelStruct.burnIn; 
isT=modelStruct.isT;

lambda=modelStruct.lambda; %shrinkage parameter
h = sqrt(1-lambda^2);
% initialise using empiricalStd or ewaStd
parts(1:burnIn,:)=randn(burnIn,N)*.1 + repmat(initVs(1:burnIn,:),[1 N]);

% unpack gpStruct
hyp = gpStruct.hyp;
% hyper-parameter matrices + df param if student-t
numParams = numel(hyp.mean)+numel(hyp.cov)+numel(hyp.lik)+isT; 
paramParts = zeros(N,numParams,T); % hypers particles
% initial hyperparam means
initHypVec = hyp2Mat(hyp,numel(hyp.mean),numel(hyp.cov),numel(hyp.lik));

learnHyp=gpStruct.learnHyp;
if learnHyp % learn using particle filters and init from prior
    paramParts(:,1:numParams-isT,1:burnIn)=randn(N,size(initHypVec,2),burnIn) ...
        + repmat(initHypVec,[N 1 burnIn]);
else % for known hyper-parameters
    paramParts(:,1:numParams-isT,1:burnIn)=repmat(initHypVec,[N 1 burnIn]);
end
if isT
    paramParts(:,end,1)=randn(N,1); %df nu-2>0 -> log(nu-2)~N()
end

% store predictive likelihoods
predLikAll = zeros(T,N); %predictive likelihood for each particle
weights=zeros(T,N);

% rapcf stats
probs = [.05 .5 .95]; %quantile probs
partQuants=zeros(T,length(probs),d);
numEffs= zeros(T,1); %tracks effective number of particles
timeUsed=zeros(T,1);

for t=burnIn+1:T-1 
    tStart=tic;
    
    prevDensity=weights(t-1,:)';
    % shrink params: gaussian diffusion
    prevParams = paramParts(:,:,t-1);
    paramEst= rapfPointEstimate(prevParams,prevDensity,lambda); 
    
    % expected log variances under gp and shrunken hypers
    partsMu = gpVolPredFun(parts(1:t-1,:),data(1:t-1,:),gpStruct,paramEst,isT);
    
    % pdfs on expected 
    mDensityMu = modelStruct.obsfunc(data(t,:),partsMu,paramEst,isT,modelStruct.transfunc);
    g=mDensityMu+prevDensity;
    
    % resample - importance weight 
    [normWgts Neff] = normalizeWeights(g);
    [pIndex]= resample(normWgts); % resample index

    % only propagate chain of IS sampled indices forward
    parts(1:t-1,:) = parts(1:t-1,pIndex);
    
    % jitter hypers
    if learnHyp
        paramsStd = repmat(std(prevParams),N,1);
        newParams = randn(size(paramEst)).*paramsStd*h + paramEst(pIndex,:); %resampled params
    else % assume fixed - don't jitter
        newParams = paramParts(:,:,t-1); % just use last
        paramParts(:,:,t) = newParams;
    end

    % debug statements
    numEffs(t,:)=Neff;
    fprintf('t=%i; numEffs=%.1f\n',t,numEffs(t));
    %inspect particle representation after resampling
    partQuants(1:t-1,:)=quantile(parts(1:t-1,:),probs,2); 
    
    % propose log variances from gp and corresponding hypers
    newParts = gpVolPredFun(parts(1:t-1,:),data(1:t-1,:),gpStruct,newParams,isT);
    
    % observation probability on new particles
    mDensity = modelStruct.obsfunc(data(t,:),newParts,newParams,isT,modelStruct.transfunc);
    weights(t,:)= mDensity-mDensityMu(pIndex,:); %adjust for proposalDensity (aligned) and save for next time-step
    
    % store particles
    parts(t,:)=newParts;
    paramParts(:,:,t)=newParams;    
    
    % predict 1- step forward
    predParts = gpVolPredFun(parts(1:t,:),data(1:t,:),gpStruct,newParams,isT);
    % prediction of t+1 at t
    predLikAll(t,:)=modelStruct.obsfunc(data(t+1,:),predParts,newParams,isT,modelStruct.transfunc);

    timeUsed(t)=toc(tStart);
end

% compute predictive likelihood given IS weights
normWeights = normalizeWeights(weights,2);
predLik = logSum(predLikAll,normWeights,0);

stats.timeUsed=timeUsed;
stats.numEffs=numEffs;
stats.partQuants=partQuants;

return;

%% debug stuff
simVs=varargin{1};  %log variances
truePredLik = log(normpdf(data,0,exp(simVs/2)));
truePredLik(1:end-1,:)=truePredLik(2:end,:); %allign with predLik indices

compLik = [truePredLik predLik];
compLik = compLik(burnIn+1:end,:);
nanmean(compLik)
plot(cumsum(compLik)./repmat(1:size(compLik,1),size(compLik,2),1)')
title('avg pred liks'); legend('truth','GP-Vol')

% particle representation
plot([simVs partQuants]); grid on;
legend('truth V', 'V 5%', 'V 50%', 'V 95%');

simHyp=varargin{2};
simHyps = hyp2Mat(hyp,numel(simHyp.mean),numel(simHyp.cov),numel(simHyp.lik));
for i=1:numParams %slightly inaccurate as not including IS weights
    paramQuants = quantile(squeeze(paramParts(:,i,:)),probs,1);
    plot([repmat(simHyps(i),[T 1]) paramQuants' ]);
    pause
end