%% demo file
% must add gp package to path
addpath(genpath('../../../gpml3.1/'))

%% 1. generate data according to gp-vol model
numStates=200;
gpFuncs.covfunc=@covSEiso; hyp.cov=[0;0];
gpFuncs.likfunc = @likGauss; sn = .1; hyp.lik = log(sn);
gpFuncs.meanfunc = {@meanSum, {@meanLinear, @meanConst}};
hyp.mean = [.99; -.1; 0]; 

obsfunc={@randn}; %normal observation model
% obsfunc={@trnd,3} %student-t observation model
[simData simVs simErrs] = gpVolGenerator(numStates,gpFuncs,hyp,obsfunc);

% visualize simulated data
figure(1);clf; plot(simData);title('data'); mean(simData)
figure(2);clf; plot(exp(simVs)); title('variance'); 

%% 2. if latent states are known can learn hypers with GP methods
y = simVs(2:end,:);
x = [simVs(1:end-1,:), simData(1:end-1,:)];
hypOut = minimize(hyp, @gp, -100, @infExact, gpStruct.meanfunc, gpStruct.covfunc, gpStruct.likfunc, x, y);

%% 3. learn gp-vol with RAPCF
% but latent variables are not observed directly and are unknown
% therefore we need MCMC methods to learn latent variables 
% as well as model hyperparameters
% we can use standard MCMC, PGAS or RAPCF
N=1000;
modelStruct.obsfunc=@gpVolObsModel;
modelStruct.isT=0;
modelStruct.inputfunc=@gpVolInput; %remove
modelStruct.transfunc=@exp;
modelStruct.burnIn=50;
modelStruct.isT=0; %1 for student-t observation model
modelStruct.lambda=.99; %shrinkage parameter

% initialize with empirical std over window
initWinSize=5;
stdWinSize=5; % tweak to make initVs and simVs similar 
initSigmas = empiricalStd(simData,initWinSize,stdWinSize);
initVs = log(initSigmas.^2); % log variance
% initVs should look similar to simVs for quick convergence
plot([initVs simVs]) 

% learn hyp on initial data
gpStruct.covfunc=@covSEiso; hyp.cov=[0;0];    
gpStruct.likfunc = @likGauss; sn = .1; hyp.lik = log(sn);
gpStruct.meanfunc = {@meanSum, {@meanLinear, @meanConst}};
hyp.mean = [.98; -.1; 0]; 
gpStruct.hyp=hyp;
gpStruct.learnHyp=1; %0 for known hypers (equivalent to running Auxiliary Chain Particle Filter)

burnIn=50;
[gpParts paramParts predLikGP] = racpf(simData,N,modelStruct,gpStruct,initVs, simVs, hyp);

% 4. pgas - code to be commented and uploaded