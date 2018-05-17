function [empStd] = empiricalStd(returnData,initWin,winSize)
[T D] = size(returnData);
if nargin<3
    winSize=T;
end
if nargin<2
    initWin=5; %some small window size
end

%initialize
empStd = nan(T,D);
initStd = nanstd(returnData);
empStd(1:initWin)=initStd;

for i=initWin+1:T
    empStd(i)=nanstd(returnData(max(1,i-winSize+1):i));
end
return

%% debug visualizations
subplot(2,1,1)
plot(empStd,'-o'); grid on; title('Empirical Std')
subplot(2,1,2)
plot(returnData,'-o'); grid on; title('Empirical Std')