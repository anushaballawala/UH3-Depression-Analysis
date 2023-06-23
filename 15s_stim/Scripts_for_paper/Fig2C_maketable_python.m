%% combine all. 

clear 
PatientID = 'DBSTRD002'; 
DBS_target = 'lVCVS'; 
roi = 'acc'; 
%% 
files = dir(sprintf('/Users/anushaallawala/Data/DBSTRD/15s_stim/data/%s_*_%s.mat',roi,DBS_target)); 
alldata = {} ;
for i = 1:numel(files)
    data{i} = load(files(i).name); 
    switch roi 
        case 'acc'
                alldata{i} = data{i}.acc ;  %%%% change roi name here. 
        case 'amy'
                alldata{i} = data{i}.amy ;  %%%% change roi name here. 
        case 'dpf'
                alldata{i} = data{i}.dpf ;  %%%% change roi name here. 
        case 'mof'
                alldata{i} = data{i}.mof ;  %%%% change roi name here. 
        case 'lof'
                alldata{i} = data{i}.lof ;  %%%% change roi name here. 
        case 'vpf'
                alldata{i} = data{i}.vpf ;  %%%% change roi name here. 
    end 
    
    all_foi_labels{i} =  data{i}.foi ; 
    all_cond_labels{i} = data{i}.condition ; 
end 

%% 

all_foi_labels = vertcat(all_foi_labels{3},all_foi_labels{6},all_foi_labels{1},...
    all_foi_labels{2},all_foi_labels{5},all_foi_labels{4}); 

all_cond_labels = vertcat(all_cond_labels{3},all_cond_labels{6},all_cond_labels{1},...
    all_cond_labels{2},all_cond_labels{5},all_cond_labels{4}); 

alldata = vertcat(alldata{3},alldata{6},alldata{1},...
    alldata{2},alldata{5},alldata{4}); 

%all_foi_labels = cat(1,all_foi_labels{:}); 
%all_cond_labels = cat(1,all_cond_labels{:}); 
%alldata = cat(1,alldata{:}); 

%% make table. 

tbl = table(alldata,'VariableNames',{'PSD'}); 

tbl{:,2} = all_foi_labels; 
tbl{:,3} = all_cond_labels; 
tbl = renamevars(tbl,["Var2","Var3"],["FOI","Condition"]); 

%highgamma
filename = sprintf('/Users/anushaallawala/Python_projects/Stim_paper/%s_%s_%s.csv',PatientID,DBS_target,roi); 
writetable(tbl,filename); 
