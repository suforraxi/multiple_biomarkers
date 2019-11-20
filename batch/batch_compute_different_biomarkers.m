%%


function batch_compute_different_biomarkers(bioName)

% setting up the path


if(strcmp(bioName,'sdDTF'))
    eeglab;
end

%% all temporal and extra temporal

cfgBatch.inDir_data       = '/home/matteo/Desktop/tle_e/converted/'; 
cfgBatch.subj_info_F      = '/home/matteo/Desktop/tle_e/info/info.tsv';
outrootFolder             = '/home/matteo/Desktop/tle_e/zscore_notch/2Dbip/';


% parameters to select
cfgBatch.cutLast          = 1;
cfgBatch.trials           = 'all';
cfgBatch.length           = 60; %seconds of new trials
cfgBatch.overlap          = 0;

cfgBatch.cutTrials        = 1;
cfgBatch.trials_ct        = 'all';
cfgBatch.length_ct        = 5; %seconds of new trials
cfgBatch.overlap_ct       = 0;

% detrend and demean
cfgBatch.deT_deM           = 1;

% notch filter
cfgBatch.notch   = 1;
cfgBatch.notchBS = [49 51];

switch bioName

    case 'ARR'
        %% ARR (Geertsema 2017)
        
        cfgBatch.montage      = 'bipolar_two_directions';
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_ARR_2Dbip.txt');
        
        cfgBatch.epiBio  = 'ARR';
        
        cfgBatch.windowL = 40; 

    case 'PAC'        
        
        cfgBatch.montage      = 'bipolar_two_directions';
  
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_PAC_2Dbip.txt');
        
        cfgBatch.epiBio  = 'PAC';
        
        cfgBatch.lb      = [4 8];
        cfgBatch.hb      = [30 80];

    case 'PLI'
        %% PLI (Stam 2007)
        cfgBatch.montage      = 'bipolar_two_directions';
        
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_PLI_2Dbip.txt');
        
        cfgBatch.epiBio  = 'PLI';
        
        cfgBatch.boi     = [30 80]; 

        compute_all_situations_bipolar2D(cfgBatch)
    case 'PLV'
        %% PLV (Mormann 2000)

        cfgBatch.montage      = 'bipolar_two_directions';
        
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_PLV_2Dbip.txt');
        
        cfgBatch.epiBio  = 'PLV';
        
        cfgBatch.boi     = [30 80]; 

    case 'H2'
        %% H2 (Kalitzin 2006)
        
        cfgBatch.montage      = 'bipolar_two_directions';

        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_H2_2Dbip.txt');
        
        cfgBatch.epiBio  = 'H2';
        
        cfgBatch.boi     = [30 80];
        cfgBatch.T       = -68:17:68;
        cfgBatch.n       = 100;
        %cfgBatch.N       = [];
        

    case 'GC'
        %% GCtime (Park 2018)
        cfgBatch.montage      = 'bipolar_two_directions';

        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_GC_2Dbip_artfree.txt');
        cfgBatch.epiBio       = 'GC';

        cfgBatch.notch         = 0;
        cfg.notchBS            = [];
        cfgBatch.momax         = 30; % max model order (AIC / BIC) 
        
    case 'sdDTF'
        %% sdDTF (Zweiphenning 2019)
       
        cfgBatch.montage      = 'bipolar_two_directions';
        
        cfgBatch.outdir_combi = fullfile(outrootFolder,'/combined/');
        cfgBatch.errorFile    = fullfile(outrootFolder,'/info/errors_sdDTF_2Dbip.txt');
        cfgBatch.epiBio       = 'sdDTF';

        cfgBatch.morder            = 30;
        cfgBatch.freqBand          = 30:80;
        cfgBatch.WindowLengthSec   = 5;
        cfgBatch.WindowStepSizeSec = 5;

        
end

compute_all_situations_bipolar2D(cfgBatch)
