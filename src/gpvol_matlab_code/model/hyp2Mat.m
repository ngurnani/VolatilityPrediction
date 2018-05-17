function [paramVec] = hyp2Mat(hyp, meanParams,covParams,likParams)
paramVec=zeros(1,meanParams+covParams+likParams);
paramVec(1:meanParams)=hyp.mean;
paramVec(meanParams+1:meanParams+covParams)=hyp.cov;
paramVec(meanParams+covParams+1:meanParams+covParams+likParams)=hyp.lik;