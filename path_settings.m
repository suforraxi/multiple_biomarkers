% 
% Setting path for all the dependencies for manuscript DOI:
% All the matlab functions and toolboxes required
% They are freely available on the web

git_root        = '/home/matteo/Desktop/git_rep/';
root_repository = fullfile(git_root,'multiple_biomarkers');
eeglab_root     = '/home/matteo/Desktop/mat_lib/eeglab2019_0/';
mvgc_root       = '/home/matteo/Desktop/mat_lib/mvgc_v1.0';

% functions to import BIDS structures for ioECoG  (https://github.com/suforraxi/ieeg_respect_bids)         
addpath(fullfile(git_root,'ieeg_respect_bids'))
addpath(fullfile(git_root,'ieeg_respect_bids','external'))
addpath(fullfile(git_root,'ieeg_respect_bids','importBIDS'))
addpath(fullfile(git_root,'ieeg_respect_bids','trc2bids'))
addpath(fullfile(git_root,'ieeg_respect_bids','micromed_utils'))

% json functions (https://github.com/fangq/jsonlab)
addpath(fullfile(git_root,'jsonlab')) 

% project functions 
addpath(fullfile(root_repository,'biomarker_computation'))
addpath(fullfile(root_repository,'biomarker_computation','biomarker_wrappers'))
addpath(fullfile(root_repository,'biomarker_computation','montage'))
addpath(fullfile(root_repository,'biomarker_result_tables'))
addpath(fullfile(root_repository,'analysis_and_figures'))
addpath(fullfile(root_repository,'batch'))
addpath(fullfile(root_repository,'read_bids'))
addpath(fullfile(root_repository,'external'))
   

% MVGC toolbox  (https://users.sussex.ac.uk/~lionelb/MVGC/)      
run(fullfile(mvgc_root,'startup.m'));

        

% eeglab (https://sccn.ucsd.edu/eeglab/index.php)

%addpath(sprintf('%s%splugins%sfirfilt2.3',eeglab_root,filesep,filesep))
%addpath(sprintf('%s%splugins%sdipfit3.2',eeglab_root,filesep,filesep))
%addpath(sprintf('%s%splugins%sBiosig3.3.0%sbiosig%st250_ArtifactPreProcessingQualityControl',eeglab_root,filesep,filesep,filesep,filesep))
%addpath(sprintf('%s%splugins%sBiosig3.3.0%sbiosig%st200_FileAccess',eeglab_root,filesep,filesep,filesep,filesep))
%addpath(sprintf('%s%splugins%sBiosig3.3.0',eeglab_root,filesep,filesep))
addpath(sprintf('%s%splugins',eeglab_root,filesep))
addpath(sprintf('%s%sfunctions%ssigprocfunc',eeglab_root,filesep,filesep))
addpath(sprintf('%s%sfunctions',eeglab_root,filesep))
addpath(sprintf('%s',eeglab_root))


% fieldtrip (https://github.com/fieldtrip)
addpath((fullfile(git_root,'fieldtrip')));

% Violin plots (https://github.com/bastibe/Violinplot-Matlab)
addpath(fullfile(git_root,'Violinplot-Matlab'));