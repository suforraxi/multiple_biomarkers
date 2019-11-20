%% Setting path for all the dependencies
%

git_root        = '/home/matteo/Desktop/git_rep/';
root_repository = fullfile(git_root,'multiple_biomarkers');
       
eeglab_root     = '/home/matteo/Desktop/mat_lib/eeglab2019_0/';
        
        
        
        
addpath(fullfile(git_root,'ieeg_respect_bids/'))
addpath(fullfile(git_root,'ieeg_respect_bids/external/'))
addpath(fullfile(git_root,'ieeg_respect_bids/importBIDS/'))
addpath(fullfile(git_root,'ieeg_respect_bids/trc2bids/'))
addpath(fullfile(git_root,'ieeg_respect_bids/micromed_utils/'))
        
addpath(fullfile(root_repository,'biomarker_wrappers/'))
addpath(fullfile(root_repository,'batch/'))
addpath(fullfile(root_repository,'episign/batch/extra_temporal/'))
addpath(fullfile(root_repository,'episign/create_montage/'))
addpath(fullfile(root_repository,'episign/bids/'))
        

        
run('/home/matteo/Desktop/mat_lib/mvgc_v1.0/startup.m');

        
addpath(fullfile(git_root,'Violinplot-Matlab/'));
