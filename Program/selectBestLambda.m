function lambda = selectBestLambda(lambdaParams, relErrMat, sumDeltaMat)
%% lambda = selectBestLambda(lambdaParams, relErrMat, sumDeltaMat)
% Select the best lambda value after cross validation.
% 1) Perform min-max scaling for relative errors and sums of corrections
% 2) Calculate the average product of scaled relative errors and sums of
%    corrections
% 3) Determine the second numerical gradient of the obtained scores in step
%    two
% 4) Find the first sign change starting from the largest lambda value and
%    return the associated lammbda value
% Input:
%   double lambdaParams:        explored lambda values
%   double relErrMat:           matrix of relative errors from
%                               cross-validation (#iterations x #lambdas)
%   double sumDeltaMat:         matrix of sums of corrections (deltas) from
%                               cross-validation (#iterations x #lambdas)
% Output:
%   double lambda:              selected optimal lambda

% min-max scaling
err_standard = (relErrMat-min(relErrMat)) ./ range(relErrMat);
log_delta_standard = (log10(sumDeltaMat)-min(log10(sumDeltaMat))) ./ range(log10(sumDeltaMat));
% average product of scaled values
lambdaScores = mean(err_standard.*log_delta_standard, 2, 'omitnan');
% numerical gradient
d1 = gradient(lambdaScores, log10(lambdaParams));
d2 = gradient(d1, log10(lambdaParams));
ifp_pos_1 = strfind(sign(d2'), [-1 1])+1;
% select best lambda
lambda = lambdaParams(ifp_pos_1(end));


end