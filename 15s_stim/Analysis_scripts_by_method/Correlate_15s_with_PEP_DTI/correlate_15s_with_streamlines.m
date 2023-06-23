%%% linear regression for # of streamlines 

%% load averaged data saved for each condition 
    %list experiment info for each channel 
    PatientID = 'DBSTRD001'; DBS_target = 'SCC'; hemi = 'l'; stim_freq = 130; 
    experiment_name = sprintf('%s%s_f%d',hemi,DBS_target,stim_freq); 

switch DBS_target 
    case 'SCC' 
        config_names = {'elec1','elec234','elec25','elec36','elec47',...
                       'elec567','elec8'}; 

    case 'VCVS' 
        config_names =  {'elec1','elec25','elec36','elec47',...
                       'elec8'};         
end 

%% load all FOI              
data_15s = load(sprintf('E:/DBSTRD/%s/Experiments/15s_stim/Processed Data/data for Logit/AllCh_allFOI_%s.mat',...
    PatientID,experiment_name));  %all FOI

% combine all FOI & theta, alpha, beta and &delta separately 

%% load streamlines data 

streamlines_data = load(sprintf('E:/DBSTRD/%s/Experiments/PEP_data/streamlines_%s%s.mat',...
    PatientID,hemi,DBS_target));

%streamlines data is already formatted to match channel labels with 15-s &
%PEP data 

%format 15s data to match indices for PEP data 
stim_15s_idx = [1:12,14:24,26:37,39:58,62:75,77:105,107:117,120:123,125:131]; 
ch_labels = data_15s.ch_labels; 
labels_15s = ch_labels(stim_15s_idx); 

%% 

for i = 1:7   
    data_15s.tbl_allch_alldata_stimpre_alpha{i, 1} = data_15s.tbl_allch_alldata_stimpre_alpha{i, 1}(stim_15s_idx);
    data_15s.tbl_allch_alldata_stimpre_beta{i, 1} = data_15s.tbl_allch_alldata_stimpre_beta{i, 1}(stim_15s_idx);
    data_15s.tbl_allch_alldata_stimpre_delta{i, 1} = data_15s.tbl_allch_alldata_stimpre_delta{i, 1}(stim_15s_idx);
    data_15s.tbl_allch_alldata_stimpre_theta{i, 1} = data_15s.tbl_allch_alldata_stimpre_theta{i, 1}(stim_15s_idx);

end 

%% create matrix that is contact configuration X channel for each FOI for 15-s stim data 

% Convert table to matrix 
alpha_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_alpha);  
theta_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_theta); 
beta_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_beta); 
delta_15s = cell2mat(data_15s.tbl_allch_alldata_stimpre_delta); 


switch DBS_target 
    case 'SCC' 
        %index out conditions that we dont have for PEP (ring configurations) 
        alpha_15s = vertcat(alpha_15s(1,:), alpha_15s(3:5,:), alpha_15s(7,:)); 
        theta_15s = vertcat(theta_15s(1,:), theta_15s(3:5,:), theta_15s(7,:)); 
        beta_15s = vertcat(beta_15s(1,:), beta_15s(3:5,:), beta_15s(7,:)); 
        delta_15s = vertcat(delta_15s(1,:), delta_15s(3:5,:), delta_15s(7,:)); 
        disp(size(alpha_15s))
    case 'VCVS'
        disp('do nothing') 
end 

%% 
vec = @(x) x(:);
X = [alpha_15s(:), theta_15s(:), beta_15s(:), delta_15s(:)];
lm = fitlm(X, log(1 + vec(streamlines_data.streamline_matrix'))); 
lm
xs = [4];
y = 2;
lm2 = fitlm(X(:,xs), X(:,y));
lm2
%% 
scatter(vec(streamlines_data.streamline_matrix'),delta_15s(:))

