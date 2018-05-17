function [simData Vs simErrs] = gpVolGenerator(numStates,gpStruct,hyp,datafunc)
simData= zeros(numStates,1); %1D 
Vs = zeros(numStates,1); %log variance
simErrs= zeros(numStates,1);

meanfunc=gpStruct.meanfunc;
covfunc=gpStruct.covfunc;
likfunc=gpStruct.likfunc;

% initialise data
SigInit=1; % var_0=1;
Vs(1)=log(SigInit); 
simData(1)=feval(datafunc{:},1,1).*exp(Vs(1)).^.5; %adjust by std

%GP-Vol model: v_t = f(v_{t-1},x_{t-1})
for t=2:numStates
    input = [Vs(1:t-2), simData(1:t-2)]; % training inputs
    target = Vs(2:t-1,:); % training targets
    z = [Vs(t-1), simData(t-1)]; % test input for next prediction
    [m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc,input,target, z);
    
    Vs(t) = m + randn(1)*(s2^.5); % should have process noise
    
    simData(t) = feval(datafunc{:},1,1)*exp(Vs(t)).^.5; %generate in std space
end