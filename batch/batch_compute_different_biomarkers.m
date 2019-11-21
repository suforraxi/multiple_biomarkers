% Functions for manuscript DOI:
%
% Batch to compute different biomarkers

% input 
%
% inDir_data     -   input folder for where the raw data in BIDS is located
% 
% subj_info_F    -   table with subject descriptions,it should contain the following variables
%
%                              subjID                 - Coded name from the database (RESPXXXX)
%                              available	          - 1 if the subject can be used 0 if the subject
%                                                       cannot be used
%                              primary_path_class	  - index from 1 to 10 indicating the primary
%                                                       pathology
%                                                       1  High Grade Tumor (WHO III + IV)
%                                                       2  Low Grade Tumor (WHO I + II)
%                                                       3  MTS
%                                                       4  FCD
%                                                       5  no abnormalities
%                                                       6  cavernoma
%                                                       7  gliosis/scar
%                                                       8  AVM
%                                                       9  malformation cortical development
%                                                       10 TuberoSclerosis
% 
%                             description_sf_1y      - seizure freedom outcome one year after surgery (Engel class + medication
%                                                     level after surgery) 
%                                                     i.e. 
%                                                      1A_AED_eq   if the subject is Engel 1A class with
%                                                                  same medication after surgery
%                                                      1A_AED_stop if the subject is Engel 1A class who
%                                                                  stopped medication after surgery
%                                                      1A_AED_low  if the subject is Engel 1A class who 
%                                                                  decreased medication after surgery
%                             description_sf_longest - longest seizure freedom outcome reported after surgery  
%                                                     (Engel class + medication level after surgery) 
% 
% 
% 
%                             typeEPI                - type of epilepsy for the subject
%                                                      T Temporal
%                                                      E Extra Temporal
% 
%                             HFOstudy               - >0 if the subject was in the HFO trial 
%                                                       0 otherwise
%
%
%  bioName      - biomarker to be computed
%                    it could be one of the following string
%                       ARR   AutoRegressive model Residual  (Geertsma 2017)
%                       PAC   Phase Amplitude Coupling 
%                       PLI   Phase Lag Index  (Stam 2007)
%                       PLV   Phase Locking Value (Mormann 2000)
%                       H2    Non linear correlation coefficient (Kalitzin 2006)
%                       GC    Time-based Granger Causality (Lionel Barnett and Anil K. Seth, 2014 MVGC Toolbox)
%                       sdDTF Short-time direct Directed Transfer Function (Mullen 2014 SIFT toolbox)
%
% outrootFolder - The results are saved in the folder specified by outrootFolder see below
%                     The generic result is saved in a matlab structure named outres with the following fields
% 
%                           outres.hdr          - containing the datasetName from which the
%                                                 biomarker is computed
%                           outres.label        - cell array containing the channel names for which
%                                                 the biomarker is computed 
%                           outres.time         - cell array containing the time for each trial in
%                                                 which the original data is divided
%                           outres.fsample      - sample frequency of the data
%                           outres.sampleinfo   - sample information corresponding to the samples
%                                                 in the original recordings from which each trial
%                                                 computed
%                           outres.bio          - cell array with an entry per trial of the input
%                                                 data. Each entry corresponds another cell array
%                                                 with the computation of the biomarker per each
%                                                 channel specified by label
%                                                 (or biomaker statistic, like the strenght for bivariate/ multivariate methods)
%                           outres.extra        - struct with extra results depending on the
%                                                 biomaker (see biomarker_wrapper specific function)
%                           outres.type         - name of the biomaker computed

function batch_compute_different_biomarkers(inDir_data,subj_info_F,outrootFolder,bioName)


if(strcmp(bioName,'sdDTF'))
    eeglab;
end

% input and output folders

% input folder for where the raw data in BIDS is located
cfgBatch.inDir_data       = inDir_data; 
cfgBatch.subj_info_F      = subj_info_F;


% parameters for the selection of the last 60 seconds of each recording
% this choice was made to reduce the propofol influence
cfgBatch.cutLast          = 1; 
cfgBatch.trials           = 'all';
cfgBatch.length           = 60; %seconds of new trials
cfgBatch.overlap          = 0;

% redefinition of the last minute in smaller trials (5 secs)
% this choice as an 'harmonization' step among the following papers
%
% Geertsema et al. 2015; 
% Geertsema et al. 2017; 
% Amiri et al. 2016; 
% Varatharajah et al. 2018; 
% Cimbalnik et al. 2019; 
% Guirgis et al. 2015; 
% Mormann et al. 2000; 
% van Dellen et al. 2009; 
% Van Diessen et al. 2013; 
% Bettus et al. 2008; 
% Park and Madsen 2018; 
% Zweiphenning et al. 2019 
%
cfgBatch.cutTrials        = 1;
cfgBatch.trials_ct        = 'all';
cfgBatch.length_ct        = 5; %seconds of new trials
cfgBatch.overlap_ct       = 0;

% detrend and demean
cfgBatch.deT_deM           = 1;

% notch filter
cfgBatch.notch   = 1;
cfgBatch.notchBS = [49 51]; % line noise at 50Hz

switch bioName

    case 'ARR'
        % ARR (Geertsema 2017)
        
        cfgBatch.montage      = 'bipolar_two_directions';
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_ARR_2Dbip.txt');
        
        cfgBatch.epiBio  = 'ARR';
        
        cfgBatch.windowL = 40; % window to compute ARR

    case 'PAC'        
        % Phase Amplitude Coupling
        cfgBatch.montage      = 'bipolar_two_directions';
  
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_PAC_2Dbip.txt');
        
        cfgBatch.epiBio  = 'PAC';
        
        cfgBatch.lb      = [4 8];  % low band where to compute phase
        cfgBatch.hb      = [30 80];% high band where to compute amplitude  

    case 'PLI'
        % PLI (Stam 2007)
        cfgBatch.montage      = 'bipolar_two_directions';
        
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_PLI_2Dbip.txt');
        
        cfgBatch.epiBio  = 'PLI';
        
        cfgBatch.boi     = [30 80]; % PLI for band of interest

        compute_all_situations_bipolar2D(cfgBatch)
    case 'PLV'
        % PLV (Mormann 2000)

        cfgBatch.montage      = 'bipolar_two_directions';
        
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_PLV_2Dbip.txt');
        
        cfgBatch.epiBio  = 'PLV';
        
        cfgBatch.boi     = [30 80];  % PLV for band of interest

    case 'H2'
        % H2 (Kalitzin 2006)
        
        cfgBatch.montage      = 'bipolar_two_directions';

        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_H2_2Dbip.txt');
        
        cfgBatch.epiBio  = 'H2';
        
        cfgBatch.boi     = [30 80];   % H2 for band of interest
        cfgBatch.T       = -68:17:68; % delays where to compute H2 in samples 
        cfgBatch.n       = 100;       % number of bins see (external/h2delay.m)
      

    case 'GC'
        % GCtime (Park 2018)
        cfgBatch.montage      = 'bipolar_two_directions';

        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_GC_2Dbip_artfree.txt');
        cfgBatch.epiBio       = 'GC';

        cfgBatch.notch         = 0;
        
        cfgBatch.momax         = 30; % max model order (AIC / BIC) 
        
    case 'sdDTF'
        % sdDTF (Zweiphenning 2019)
       
        cfgBatch.montage      = 'bipolar_two_directions';
        
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_sdDTF_2Dbip.txt');
        cfgBatch.epiBio       = 'sdDTF';

        cfgBatch.morder            = 30;    % model order according to (Zweiphenning 2019)
        cfgBatch.freqBand          = 30:80; % band of interest
        cfgBatch.WindowLengthSec   = 5;     % window size for trial
        cfgBatch.WindowStepSizeSec = 5;     % no overlapping

        
end

compute_all_situations_bipolar2D(cfgBatch)
