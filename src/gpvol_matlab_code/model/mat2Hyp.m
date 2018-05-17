function [newHyp] = mat2Hyp(paramVec,meanParams,covParams,likParams)
[a p] =size(paramVec);
if p~=meanParams+covParams+likParams
    error('incorrect hyp dimensions')
end
newHyp.mean=paramVec(1:meanParams)'; %col vectors
newHyp.cov=paramVec(meanParams+1:meanParams+covParams)'; %col vectors
newHyp.lik=paramVec(meanParams+covParams+1:meanParams+covParams+likParams)';