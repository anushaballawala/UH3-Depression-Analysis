%% Make network map 

addpath(genpath('/Users/anushaallawala/Documents/MATLAB/Github/Packages/CircularGraph/'))

label_tbl = ROI_labels.summary_ch_info;
gray_mt_idx = find(contains(label_tbl.Matter,'Grey')==1); 
label_tbl_new = ROI_labels.summary_ch_info(gray_mt_idx,:);
bolt_idx = find(contains(label_tbl_new.ROI,'Bolt')==1);
label_tbl_new(bolt_idx,:) = [];

%% 


DBStarget = 'rVCVS';

switch DBStarget
    case 'lSCC'
        X = data_lSCC_diff(gray_mt_idx,gray_mt_idx,:); 
    case 'rSCC' 
        X = data_rSCC_diff(gray_mt_idx,gray_mt_idx,:); 
    case 'lVCVS'
        X = data_lVCVS_diff(gray_mt_idx,gray_mt_idx,:); 
    case 'rVCVS'
        X = data_rVCVS_diff(gray_mt_idx,gray_mt_idx,:); 
end 

%% 


        %bolt_idx = [25,38,46,73,115,116];
        X(bolt_idx,:,:) = [];
        X(:,bolt_idx,:) = [];

        %bolt_idx = [40,99];


        %bolt_idx = [27,46]; 


%% Choose threshold (maybe customize this since some values are negative) 
% thresh = 0.1;
% X(X >  thresh) = 1;
% X(X <= thresh) = 0;
% X(isnan(X)) = 0; 
%% If you want to show line thickness, customize numbers. 

thresh = 0.05; 
X(X <= thresh) = 0;
X(X > thresh) = 1;

% The values of X are initially in any arbitrary range.
% When x < -1 --> set to 0
% When -1 <= x <= -0.05 --> set to 1
% When x > -0.05 --> set to 0
% good_idx = (-1 <= X) & (X <= -thresh);
% X(good_idx) = 1;
% X(~good_idx) = 0;

%decreased coherence 

% thresh = 0.05 
% X(X <= thresh) = 0; 
% X(-1 > X > thresh) = 1; 
%% Create figure 



myLabel = label_tbl_new.ROI;
%myLabel = label_tbl.Label;


%myColorMap = lines(length(X));

figure;
circularGraph(X(:,:,1),'Label',myLabel);
filename = sprintf('/Users/anushaallawala/Desktop/Manuscripts/DBSTarget_DiffModulation_2022/Figures/Coherence/%s_%s_delta_pos_coh_%s.svg',PatientID, DBStarget,date); 
%saveas(gcf,filename);

% %% Show decreased coherence 
% thresh = -0.08;
% X(X >  thresh) = 0;
% X(X <= thresh) = 1;
% myLabel = label_tbl.ROI;
% figure;
% circularGraph(X(:,:,1),'Label',myLabel);
% title(DBStarget) 

%%
% Click on a node to make the connections that emanate from it more visible
% or less visible. Click on the 'Show All' button to make all nodes and
% their connections visible. Click on the 'Hide All' button to make all
% nodes and their connections less visible.
