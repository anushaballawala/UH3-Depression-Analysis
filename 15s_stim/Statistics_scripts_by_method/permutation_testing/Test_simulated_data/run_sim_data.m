%% Simulate sample data to validate code on with a known probability distribution 
% fake data created to run statistical pipeline through it  to check 
%  pattern of p-values is as expected. 

%% load simulated data 


load('sim_data_15s_permtest_111721.mat'); 

sim_num_ch = size(sim_data,2);
sim_num_trials = size(sim_data,1); 
sim_num_freqbands = size(sim_data,3);
%% Run permutations on sim data 

num_perm = 4; 
[t_stat] = permutations(sim_data,labels,'Ch',num_perm,@Tfunc);

%% get p value 
M = num_perm+1; 
p_value = (squeeze(mean(t_stat>=t_stat(1,:,:))))./num_perm; 

%% 

for i = 1:sim_num_ch
    figure()
    histogram(sim_data(1:5,i,1),'BinWidth',0.1); 
    hold on 
    histogram(sim_data(6:10,i,1),'BinWidth',0.1); 
    title(sprintf('Channel %d',i)); 
    figurename = sprintf('trial_distribution_sim_data_%s_Ch%d.png',date,i); 
    saveas(gcf,figurename); 
end 

%% 

histogram 
