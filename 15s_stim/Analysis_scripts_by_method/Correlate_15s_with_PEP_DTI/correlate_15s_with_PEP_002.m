%% load averaged data saved for each condition 
    %list experiment info for each channel 
    PatientID = 'DBSTRD002'; DBS_target = 'VCVS'; hemi = 'l'; stim_freq = 130; 
    experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 
    config_names = {'elec1','elec234','elec25','elec36','elec47',...
                       'elec567','elec8'}; 
               
%% load all FOI              
data_15s = load(sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/data for Logit/AllCh_allFOI_%s.mat',...
    PatientID,experiment_name));  %all FOI

% combine all FOI & theta, alpha, beta and &delta separately 

%% load all PEP matrices 

files_PEP = dir(sprintf('E:/DBSTRD/%s/Experiments/PEP_data/Threshold/Threshold_PEP-%s%s-*.mat',...
    PatientID,hemi,DBS_target)); 

for i = 1:numel(files_PEP)
    PEP_data(i) = load(files_PEP(i).name); 
end 
% concatenate& combine  all PEP current configurations 

%% index to match channel count and labels for 15-s stim data 

% load MontageInfo 
load(sprintf('E:/DBSTRD/%s/Experiments/PEP_data/MontageInfo.mat',PatientID)); 
% load channel labels for 15s stim data 
ch_labels = data_15s.ch_labels; 
 
%% Index data to match up channel labels 

switch PatientID 
    case 'DBSTRD001' 
        stim_15s_idx = [1:12,14:24,26:37,39:58,62:75,77:105,107:117,120:123,125:131]; 
        stim_PEP_idx = [2:13,15:25,27:36,38:39,41:45,47:49,51:62,64:77,79:93,95:108,110:120,122:132]; 
    case 'DBSTRD002' 
        stim_15s_idx = [1:11,13:41,43:57,59,62:65,67:101,103:127,129]; 
        stim_PEP_idx = [2:14,16:30,32:42,44:47,49,51:54,56:63,67:80,84:97,99:109,111:123,125:137]; 
end 

labels_15s = ch_labels(stim_15s_idx); 
labels_PEP = MontageInfo.sEEG.Contacts.Labels(stim_PEP_idx); 

for i = 1:numel(labels_PEP) 
   if ~startsWith(labels_15s(i),labels_PEP{i})
       disp(labels_15s(i))
       disp(i)
       disp(labels_PEP{i})
       disp('error')
   end    
end 


if numel(stim_15s_idx) == numel(stim_PEP_idx) 
    disp('ALIGNED') 
elseif numel(stim_15s_idx) ~= numel(stim_PEP_idx)
    disp('ERROR - check number of channels and labels')
end 
%% 

for i = 1:5   
    data_15s.tbl_allch_alldata_stimpre_alpha{i, 1} = data_15s.tbl_allch_alldata_stimpre_alpha{i, 1}(stim_15s_idx);
    data_15s.tbl_allch_alldata_stimpre_beta{i, 1} = data_15s.tbl_allch_alldata_stimpre_beta{i, 1}(stim_15s_idx);
    data_15s.tbl_allch_alldata_stimpre_delta{i, 1} = data_15s.tbl_allch_alldata_stimpre_delta{i, 1}(stim_15s_idx);
    data_15s.tbl_allch_alldata_stimpre_theta{i, 1} = data_15s.tbl_allch_alldata_stimpre_theta{i, 1}(stim_15s_idx);

end 


for i = 1:7 
    PEP_data(i).PEPOccurrence = PEP_data(i).PEPOccurrence(stim_PEP_idx);  
end 

%% create matrix that is contact configuration X channel for each FOI for 15-s stim data 

% Convert table to matrix 
alpha_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_alpha);  
theta_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_theta); 
beta_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_beta); 
delta_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_delta); 


% switch DBS_target 
%     case 'SCC' 
%         %index out conditions that we dont have for PEP (ring configurations) 
%         alpha_15s = vertcat(alpha_15s(1,:), alpha_15s(3:5,:), alpha_15s(7,:)); 
%         theta_15s = vertcat(theta_15s(1,:), theta_15s(3:5,:), theta_15s(7,:)); 
%         beta_15s = vertcat(beta_15s(1,:), beta_15s(3:5,:), beta_15s(7,:)); 
%         delta_15s = vertcat(delta_15s(1,:), delta_15s(3:5,:), delta_15s(7,:)); 
%         disp(size(alpha_15s))
%     case 'VCVS'
%         disp('do nothing') 
% end 
%% Flatten matrices 

theta_15s_flattened = theta_15s(:); 
alpha_15s_flattened = alpha_15s(:); 
beta_15s_flattened = beta_15s(:); 
delta_15s_flattened = delta_15s(:); 

% flatten and add on other FOIs. 
allFOI_15s_flattened = [alpha_15s(:), beta_15s(:), delta_15s(:),  theta_15s(:)] ; 
 

%% create matrix of PEP values and flatten 

switch DBS_target 
    case 'VCVS' 
        alldata_PEP = [PEP_data(1).PEPOccurrence, PEP_data(3).PEPOccurrence ,...
    PEP_data(4).PEPOccurrence, PEP_data(5).PEPOccurrence,... 
    PEP_data(7).PEPOccurrence]'; 
    case 'SCC' 
            alldata_PEP = [PEP_data(1).PEPOccurrence, PEP_data(2).PEPOccurrence ,...
        PEP_data(3).PEPOccurrence, PEP_data(4).PEPOccurrence,... 
        PEP_data(5).PEPOccurrence,PEP_data(6).PEPOccurrence,PEP_data(7).PEPOccurrence]'; 
end 

alldata_PEP_flattened = alldata_PEP(:);


%% Perform logistic regression. 

neural_feature = input('enter neural feature of interest ','s'); 
switch neural_feature
    case 'theta' 
        X = theta_15s_flattened; 
    case 'alpha'
        X = alpha_15s_flattened; 
    case 'beta'
        X = beta_15s_flattened;
    case 'delta' 
        X = delta_15s_flattened; 
    case 'all features'
        X = allFOI_15s_flattened; % matrix that is flattened ch x num_features. ?* 
end 
cutoff_value = 0.25; 
y = alldata_PEP_flattened; %logical values of 1 or 0 for PEP 


% Fit a logistic regression to get the model coefficients.
b = glmfit(X,y,'binomial','link','logit');

yhat = glmval(b,X,'logit');

yhat_binary = yhat > cutoff_value;

% Get ROC curve 

mdl = fitglm(X,y,'Distribution','binomial','Link','logit'); 
scores = mdl.Fitted.Probability; 
[X1,Y1,T,AUC] = perfcurve(alldata_PEP_flattened,scores,'true'); 


% plot data 
f = figure(); 
f.Position = [100 100 1500 450]; 
subplot(131)
confusionchart(y,yhat_binary);
title(sprintf('%s',neural_feature)); 
% 
subplot(132) 
scatter(X,y)
title(sprintf('Scatter plot of feature vs PEP - %s',neural_feature));  

subplot(133)
plot(X1,Y1)
hold on 
hline = refline([1 0]); 
hline.Color = 'r'; 
xlabel('False Positive Rate')
ylabel('True Positive Rate') 
dim = [0.695338923829489 0.847500001266599 0.0496156521188389 0.0674999987334013]; 
str = sprintf('%s',string(AUC)); 
annotation('textbox',dim,'String',str,'FitBoxToText','on'); 

% 
% figure()
% scatter(double(y), yhat, 'ro', 'markerfacealpha',0.3, 'markeredgealpha',0.1, 'markerfacecolor','r');
% axis([-0.25, 1.25, -0.25, 1.25]);
% title(sprintf('%s',neural_feature));








