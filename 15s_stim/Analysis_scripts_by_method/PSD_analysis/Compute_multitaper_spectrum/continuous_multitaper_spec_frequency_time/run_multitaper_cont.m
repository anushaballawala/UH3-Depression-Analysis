%% run multitaper spectrum over time on time series data 

addpath(genpath('/Users/anushaallawala/Data/Baseline_Data'))

%bl_data = load('/Users/anushaallawala/Data/Baseline_Data/001/referencedsig_BaselineFix_run-07_blk-01.mat'); 
%data = bl_data.all_ref_signals; 

PatientID = 'DBSTRD001';

bl_data = load('/Users/anushaallawala/Data/Baseline_Data/001/referencedsig_BaselineFix_run-07_blk-01'); 
data = bl_data.all_ref_signals'  ; 
num_ch = size(data,1); 
%% Run chronux multitaper spectrum continuous fxn 

addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/chronux_2_12')); 

% add params for multitaper 
params.Fs = 1000; % 1KHz 
params.fpass = [1 150];
params.err = [1 0.05];
params.trialave = 0;
params.tapers = [4 7]; 

window = 0.5; % 0.5s 
winstep = 0.05; % 500 ms 
movingwin = [window winstep] ; 
S = []; t = []; 
for i = 1:num_ch
 [S(:,:,i),t,f,~] = mtspecgramc(data(:,i),movingwin,params); 
end 

%% 
plot_matrix(S(:,:,1),t,f);
xlabel([]); % plot spectrogram
%caxis([8 28]); colorbar;


