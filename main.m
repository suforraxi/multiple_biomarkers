% main script manuscript DOI:



% setting all funciton dependencies

path_settings;
% biomarker names to compute
bioNames = {'ARR','PAC','PLV','PLI','H2','GC','sdDTF'};
%bioNames = {'PAC','sdDTF'};

% flag to decide if compute biomarkers or just plot the results in case the
% data from the biomarkers was already computed
compute_bio  = 1;
save_tbl_fig = 0;

% compute and save all the biomarkers

% input folder where the raw data in BIDS is stored 
inDir_data       = '/home/matteo/Desktop/tle_e/converted/';
%inDir_data       = '/home/matteo/Desktop/analysis_multiple_biomarkers/converted/'; 


% table with information related to subjects 
% (see batch_compute_different biomakers help for a description of the table 
%subj_info_F      = '/home/matteo/Desktop/analysis_multiple_biomarkers/info/info.tsv';
subj_info_F      = '/home/matteo/Desktop/replicate_analysis/info/info.tsv';


% output folder where the results from the computation of the biomarker is
% stored (see batch_compute_different biomakers help for a description of the how the struct that is saved)

%outFolder        = '/home/matteo/Desktop/tle_e/zscore_notch/2Dbip/combined/';
%outFolder        = '/home/matteo/Desktop/analysis_multiple_biomarkers/biomarker_res/';
outFolder        = '/home/matteo/Desktop/replicate_analysis/';
% output folder used to save biomarker summary tables (see create_summary_table.m for the layout of the table) 
% and figures 
%root_outFolder   = '/home/matteo/Desktop/analysis_multiple_biomarkers/';
root_outFolder   = '/home/matteo/Desktop/replicate_analysis/';

% compute biomarkers
if(compute_bio)
    parfor i = 1 : numel(bioNames)
    %for i = 1 : numel(bioNames)

       batch_compute_different_biomarkers(inDir_data,subj_info_F,outFolder,bioNames{i});
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
    cfg.rootSummaryFolder = fullfile(root_outFolder ,'summary_tables'); 

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
                         };

    % regular expression defining the epilepsy type (all/ T = temporal, E = extra-temporal ) 
    cfg.typeEPI = {'\w*','T','E'};  

    % seizure outcome to consider (i.e. after one year/longer)
    cfg.sf_class   = {'description_sf_1y'};

    % biomarker of interest 
    cfg.bioNames = bioNames;


    % comparison between maximum distribution per subject
    cfg.bioMarker2plotMaxDistribution = [1 2 3 4 5 6 7 ];  
    cfg.alpha_level                   = 0.01            ;

    % selection of a specific pathology group (i.e. MST, FCD etc they are coded
    % in information table) or all subjects
    cfg.path_label                 = {'all'};
    cfg.path_group                 = {[0]};
    % pathology group to consider according to what group is investigated 
    cfg.path_idx_of_interest       = 1;
    % regular expression to define the seizure outcome groups where to test the global
    % threshold (pre-resection recordings in improved patients (1A|1B Engel))
    cfg.seizOut2try                = {'1a_AED_stop\w*','1(a|b)\w*'};
    cfg.seizOut_idx                = [1 2];

    %regular expression to define improved patients for the group analysis
    %pooling the channels together
    cfg.seizOut_poolingCH          = {'1(a|b)\w*'};
    
    % remove channels labeled with N-N (not present)
    cfg.removeNNchannels           = 1;

    % function to plot and save the figures
    compute_results_and_save_figure(cfg)
end