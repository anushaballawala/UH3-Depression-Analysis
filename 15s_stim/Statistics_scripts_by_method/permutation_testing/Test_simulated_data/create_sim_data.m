%% Simulate sample data to validate code on with a known probability distribution 

sim_num_ch = 5;
sim_num_trials = 10; 
sim_num_freqbands = 4;
sim_data = zeros(sim_num_trials,sim_num_ch,sim_num_freqbands); 


%% Condition 1 Stim OFF 
% Channel 1 (Big difference between 2 conds)
sim_data(1,1,:) = [5, 10, 3, 2]; 
sim_data(2,1,:) = [5.1, 10.1, 3, 2]; 
sim_data(3,1,:) = [6, 12, 4, 1.3]; 
sim_data(4,1,:) =  [7, 13, 4.1, 1.4]; 
sim_data(5,1,:) =  [6.5, 12.5 4.5, 2]; 

%Channel 2 (no diff, lots of variability b/w trials)
sim_data(1,2,:) = [1,1,1,1]; 
sim_data(2,2,:) = [10, 10, 10, 8]; 
sim_data(3,2,:) = [5, 5, 5, 5]; 
sim_data(4,2,:) = [2, 2, 2, 2]; 
sim_data(5,2,:) = [10, 10, 10, 10]; 

%Channel 3 (no diff)
sim_data(1,3,:) = [0.5, 0.5, 0.5, 0.5]; 
sim_data(2,3,:) = [0.4, 0.4, 0.4, 0.4]; 
sim_data(3,3,:) = [0.5, 0.5, 0.5, 0.5]; 
sim_data(4,3,:) = [0.4, 0.4, 0.4, 0.4]; 
sim_data(5,3,:) = [4, 0.45, 0.45, 0.45]; 

%Channel 4 (big difference)
sim_data(1,4,:) = [6, 12, 2, 2]; 
sim_data(2,4,:) = [7, 13, 3, 3]; 
sim_data(3,4,:) = [6.2, 12.5, 2.5, 2.6]; 
sim_data(4,4,:) = [7.3, 13.2, 3.2, 3.2]; 
sim_data(5,4,:) = [6, 12, 2, 2];

%Channel 5 (small difference) 
sim_data(1,5,:) = [3, 7, 1, 1]; 
sim_data(2,5,:) = [5, 7.6, 1.5, 1.5]; 
sim_data(3,5,:) = [4, 7, 1.5, 1.5]; 
sim_data(4,5,:) = [3.5, 7.5, 1.2, 1.2]; 
sim_data(5,5,:) = [4.5, 7.5, 1.7, 1.2]; 

%% Condition 2 (Stim ON) 
% Channel 1
sim_data(6,1,:) = [10,20,6,1]; 
sim_data(7,1,:) = [15,21,5,1]; 
sim_data(8,1,:) = [13,22,5.5,1.5]; 
sim_data(9,1,:) = [14,22.9,5.7,1.6]; 
sim_data(10,1,:) = [14.4,22.4,5.4,1.6]; 

%Channel 2
sim_data(6,2,:) = [0.5,0.7,0.8,0.85]; 
sim_data(7,2,:) = [0.6,0.75,0.85,0.9]; 
sim_data(8,2,:) = [0.55,0.8,0.75,0.9]; 
sim_data(9,2,:) = [0.53,0.7,0.7,0.88]; 
sim_data(10,2,:) = [0.48,0.65,0.65,0.9];

%Channel 3
sim_data(6,3,:) = [0.5,0.6,0.6,0.6]; 
sim_data(7,3,:) = [0.48,0.65,0.65,0.9];
sim_data(8,3,:) = [0.53,0.7,0.7,0.88]; 
sim_data(9,3,:) = [0.6,0.75,0.85,0.9]; 
sim_data(10,3,:) = [0.5,0.7,0.8,0.85]; 

%Channel 4
sim_data(6,4,:) = [24,20,0.5,1.5]; 
sim_data(7,4,:) = [21,19,0.6,2];
sim_data(8,4,:) = [23,18,0.4,1.6];
sim_data(9,4,:) = [22,16,0.45,2.3];
sim_data(10,4,:) = [23.5, 17.5, 0.55, 3]; 

%Channel 5
sim_data(6,5,:) = [7,10,3,3]; 
sim_data(7,5,:) = [7.5,10.5,3.5,3.5]; 
sim_data(8,5,:) = [8,11,4,4];
sim_data(9,5,:) = [9,10.5,2.8,2.9];
sim_data(10,5,:) = [9,11,3,3]; 

%% Create stim labels 


sim_labels = zeros(10,1); 

sim_labels(1:5) = 0; 
sim_labels(6:10) = 1; 




