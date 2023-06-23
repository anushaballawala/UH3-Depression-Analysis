function [pvals,acceptPW,acceptSIM] = exchangeHT(X,alpha,sflag)
% function [pvals,acceptPW,acceptSIM] = exchangeHT(X,alpha,sflag=false)
%
% acceptPW and acceptSIM do take alpha values into account (you can put in
% t-stat and get out which t-stats that are significant) 

%  computes hypothesis tests of exchangeability
%
%  pvals(j,k) is a test that X(j,:) are exchangeable
%
%  k determines the type of test:
%  k == 1: upper one-sided pointwise
%  k == 2: lower one-sided pointwise
%  k == 3: two-sided pointwise
%  k == 4: upper one-sided simultaneous
%  k == 5: lower one-sided simultaneous
%  k == 6: two-sided simultaneous
%
%  [acceptPW(j,1) , acceptPW(j,2)] is a pointwise (1-alpha)*100% acceptance
%   region for X(j,1).  X(j,1) is in this region if and only if
%   pvals(j,3) > alpha. %acceptance regions where you're not significant
%   (like t - values from 0-3) 

% will be for each permutation distribution for middle 95% or lower 95% if
% its lower test 
%
%  [acceptSIM(j,1) , acceptSIM(j,2)] is a simultaneous (1-alpha)*100%
%   acceptance region for X(j,1).  X(j,1) is in this region if and only if
%   pvals(j,6) > alpha.

%accept SIM will give you values that are significant (constant for every
%single brain region and frequncy) 
%
% sflag = true first self-standardizes each row by substracting its sample
%  mean and dividing by its sample standard deviation
% the sample mean and standard deviation are computed by leaving out the
%  maximum and minimum value

if nargin < 2, alpha = []; end
if nargin < 3, sflag = false; end

if nargout > 1 && isempty(alpha)
    error('alpha must be given to create acceptance regions')
end

if ~isempty(alpha) && (numel(alpha) ~= 1 || alpha < 0 || alpha > 1)
    error('alpha must be a scalar between 0 and 1')
end

if nargout > 1 && (alpha >= .5)
    warning(['computing ' num2str((1-alpha)*100) '\% acceptance regions'])
end

% sizing
[M,Nplus1] = size(X);

% preprocessing
if sflag %standardization, set to false) 
    mu = (sum(X,2)-max(X,[],2)-min(X,[],2))/(Nplus1-2);
    sigma = sqrt(max(((sum(X.^2,2)-max(X,[],2).^2-min(X,[],2).^2) - (Nplus1-2).*mu.^2)./(Nplus1-3),eps(X(:,1))));
    X = (X - repmat(mu,1,Nplus1))./repmat(sigma,1,Nplus1);
end

X1 = repmat(X(:,1),1,Nplus1); % %takes original from permutation matrix 

pvals = zeros(M,6); 

% compute the unadjusted (pointwise) p-values
pvals(:,1) = mean(X >= X1,2);   % upper one-sided
pvals(:,2) = mean(X <= X1,2);   % lower one-sided
pvals(:,3) = min([ones(M,1), 2*pvals(:,1), 2*pvals(:,2)],[],2);  % two-sided

% compute the adjusted (simultaneous) p-values
Xmax = max(X,[],1);
Xmin = min(X,[],1);
pvals(:,4) = mean(repmat(Xmax,M,1) >= X1,2); % computes maximum over all brain regions. 1000 different maximums one for each brain region; 
%compares first X to the maximum, and finds out how many times its bigger
%than tha. If maximum is bigger only a few times, than get small p value 
pvals(:,5) = mean(repmat(Xmin,M,1) <= X1,2);
pvals(:,6) = min([ones(M,1), 2*pvals(:,4), 2*pvals(:,5)],[],2);

% compute the pointwise 2-sided acceptance regions
if nargout > 1

    alpha = floor(alpha*Nplus1/2);
    alphaN = Nplus1-alpha+1;

    if alpha == 0
        acceptPW = repmat([-inf inf],M,1);
    else
        Xs = sort(X,2);
        acceptPW = Xs(:,[alpha alphaN]);
        if sflag
            acceptPW = acceptPW.*repmat(sigma,1,2)+repmat(mu,1,2);
        end
    end

end

% compute the simultaneous 2-sided acceptance regions
if nargout > 2

    if alpha == 0
        acceptSIM = repmat([-inf inf],M,1);
    else
        Xmax = sort(Xmax);
        Xmin = sort(Xmin);
        acceptSIM = repmat([Xmin(alpha) Xmax(alphaN)],M,1);
        if sflag
            acceptSIM = acceptSIM.*repmat(sigma,1,2)+repmat(mu,1,2);
        end

    end

end
