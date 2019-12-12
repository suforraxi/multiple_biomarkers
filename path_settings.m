% 
% Setting path for all the dependencies
% All the matlab functions and toolboxes required
% They are freely available on the web

%     Copyright (C) 2019 Matteo Demuru
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.



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

addpath(sprintf('%s%splugins',eeglab_root,filesep))
addpath(sprintf('%s%sfunctions%ssigprocfunc',eeglab_root,filesep,filesep))
addpath(sprintf('%s%sfunctions',eeglab_root,filesep))
addpath(sprintf('%s',eeglab_root))


% fieldtrip (https://github.com/fieldtrip)
addpath((fullfile(git_root,'fieldtrip')));

% Violin plots (https://github.com/bastibe/Violinplot-Matlab)
addpath(fullfile(git_root,'Violinplot-Matlab'));