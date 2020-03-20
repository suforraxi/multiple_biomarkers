% main script


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


% setting all funciton dependencies
path_settings;

% biomarker names to compute
%bioNames = {'ARR','PAC','PLV','PLI','H2','GC','sdDTF'};
bioNames = {'PLV','PLI','H2'};
bioNames = {'PAC'};
bioNames = {'PAC_theta_HFA','PAC_theta_gamma'};
bioNames = {'ARR','PAC','PLV','PLI','H2','GC','sdDTF'};

% flag to decide if compute biomarkers or just plot the results in case the
% data from the biomarkers was already computed
compute_bio  = 0;
save_tbl_fig = 1;

% compute and save all the biomarkers

% input folder where the raw data in BIDS is stored 
inDir_data       = '/home/matteo/Desktop/tle_e/converted/';



% table with information related to subjects 
% (see batch_compute_different biomakers help for a description of the table 
subj_info_F       = '/home/matteo/Desktop/mb_reviewers/info/info.tsv';%'/home/matteo/Desktop/rep_analysis/info/info.tsv';

% output folder where the results from the computation of the biomarker is
% stored (see batch_compute_different biomakers help for a description of the how the struct that is saved)
%bandDir          = '13_30';
outFolder        = fullfile('/home/matteo/Desktop/mb_reviewers/combined/combined_PAC/');%'/home/matteo/Desktop/rep_analysis/combined/';



% output folder used to save biomarker summary tables (see create_summary_table.m for the layout of the table) 
% and figures 
%
root_outFolder   = '/home/matteo/Desktop/mb_reviewers/';%'/home/matteo/Desktop/rep_analysis/';

bi_boi   = {[1 4],[4 8],[8 13],[13 30]};



% choose PAC only once;

biomarker = repmat(bioNames,1,4);
bi_boi    = repmat(bi_boi,1,3);

%biomarker = ['PAC' , biomarker];
%bi_boi    = [[0 0],bi_boi];

N = numel(biomarker);

% compute biomarkers
if(compute_bio)
    parfor i = 1 : N
        
        batch_compute_different_biomarkers(inDir_data,subj_info_F,outFolder,biomarker{i},bi_boi{i});
        
    end    
end
 


% organize results in tables, plot and save figures
if(save_tbl_fig)


    cfg = [];

    % subject information table 
    cfg.info_F           = subj_info_F;

    % root folder name of raw data in BIDS format
    cfg.bidsFolder       = inDir_data;


    % Result folders and files

    % Result folder where the computed biomarkers are saved 
    cfg.rootInResFolder   = outFolder;
    % Folder where to save the summary tables
    cfg.rootSummaryFolder = fullfile(root_outFolder ,'summary_tables','new_pool'); 

    % out filename for the comparison of distributions pre and post pooling all
    % channel together
    cfg.poolingChannelFile = fullfile(root_outFolder ,'figures','pooling_channels') ;

    % Folder for the results of the comparison between distribution of the
    % maximum value per subject
    cfg.outFolderMaxComparison   = strcat(fullfile(root_outFolder ,'figures'),filesep) ;

    % regular expression to define the group of patients used to compute the global
    % threshold (cured patients) 
    %cfg.normTempRegExp = '1a_AED_stop\w*';

    % subjects to remove from the analysis
    cfg.subj2rem       = {                                                        ...                                              
                          'RESP0448','RESP0480','RESP0482','RESP0484','RESP0497', ... % hfo trial
                          'RESP0500','RESP0519','RESP0537','RESP0542','RESP0556', ... % hfo trial
                          'RESP0566','RESP0572','RESP0604','RESP0631','RESP0636', ... % hfo trial
                          'RESP0640','RESP0644', ...                                  % hfo trial 
                          'RESP0353','RESP0623'  ...                                  % no last post recording
                          'RESP0218','RESP0301'  ...                                  % frequency sample not 2048
                          'RESP0118' ...                                              % no resected channels
                         };

    % regular expression defining the epilepsy type (all/ T = temporal, E = extra-temporal ) 
    cfg.typeEPI       = {'\w*','T','E'};  
    cfg.typeEPI_label = {'Joint_T_E','Temporal','Extra_Temporal'};
    % seizure outcome to consider (i.e. after one year/longer)
    cfg.sf_class   = {'description_sf_1y'};

    % biomarker of interest 
    cfg.bioNames = bioNames;


    % comparison between maximum distribution per subject
    cfg.bioMarker2plotMaxDistribution = [1 2 4 5 6 ];  % plot only biomarkers that are significant in the group analysis pooling channels 
   % cfg.bioMarker2plotMaxDistribution = [1 2 3 ];  % plot only biomarkers that are significant in the group analysis pooling channels 
   % cfg.bioMarker2plotMaxDistribution = [1 2];  % plot only biomarkers that are significant in the group analysis pooling channels 
   
    cfg.alpha_level                   = 0.01; % alpha level for the Kolmogorov-Smirnov test

    % selection of a specific pathology group (i.e. MST, FCD etc they are coded
    % in information table) or all subjects (see primary pathology class in create_summary_table_main)
   
    cfg.path_group                    = {[0]};
    cfg.path_group_label              = {'all_primary_pathologies'};
    % pathology group to consider according to what group is investigated 
    cfg.path_idx_of_interest          = 1;
    
    % regular expression to define the seizure outcome groups where to test the global
    % threshold (pre-resection recordings in improved patients (1A|1B Engel))
    cfg.seizOut2try                   = {'1a_AED_stop\w*','1(a|b)\w*'};
    cfg.seizOut_idx                   = [1 2];

    %regular expression to define improved patients for the group analysis
    %pooling the channels together
    cfg.seizOut_poolingCH          = {'1(a|b)\w*'};
    
    % remove channels labeled with N-N (not present)
    cfg.removeNNchannels           = 1;

    % function to plot and save the figures
    compute_results_and_save_figure(cfg)
end