% Load flower data.
load fisheriris;

% skip the first 50 sasmples just so that we have 50 of each class
first_sample = 1;
% Take some subset of the data, all features for those data points
X = meas(first_sample:end,:);
% Create the label (0 vs 1)
y = strcmp('versicolor',species(first_sample:end)); % Things that are "versicolor" or not

% Fit a logistic regression to get the model coefficients.
b = glmfit(X,y,'binomial','link','logit');

yhat = glmval(b,X,'logit');

yhat_binary = yhat > 0.5;
% cutoff, i.e. 0.5 value - probabilities are funky - if number of 1's and
% 0's are not equal, this cutoff value would need to change 
%prior probability of how many 1's are there in entire sample set, and then
%probability based on feature can be adjusted 
% if the probability we gave it was 1 (because every dataset was a 1) then 

%class balance is important 
confusionchart(y,yhat_binary);

figure();
scatter(double(y), yhat, 'ro', 'markerfacealpha',0.3, 'markeredgealpha',0.1, 'markerfacecolor','r');
axis([-0.25, 1.25, -0.25, 1.25]);