    
%% create tables from outres 

% subject information table 
cfg.info_F           = '/home/matteo/Desktop/tle_e/info/info.tsv';

% folder name of subject raw data in BIDS format
cfg.bidsFolder       = '/home/matteo/Desktop/tle_e/converted/';

% remove channels not present (not anymore necessary FIX IT!!)
cfg.removeNNchannels = 1;

% expression to evaluate 'average across the epochs '
% the abs function is used because of PAC, the complex number was saved
cfg.evalBioStr       = 'nanmean(abs(cell2mat(outres.bio(~idx_art_trial))),2)';

% Result folder
rootInResFolder   = '/home/matteo/Desktop/tle_e/zscore_notch/2Dbip/combined/';
rootSummaryFolder = '/home/matteo/Desktop/tle_e/zscore_notch/2Dbip/summary_table/';


inFolder = dir(rootInResFolder);
inFolder = inFolder(~ismember({inFolder.name},{'.','..'}));


for d = 1 : numel(inFolder)

    folderName     = fullfile(inFolder(d).folder,inFolder(d).name);    
    cfg.resFolder  = folderName;

    cfg.fName2save = fullfile(rootSummaryFolder,strcat('summary_tbl_',inFolder(d).name,'.mat'));
    
    %create_summmary_table_main(cfg)

end

%% max tables

% outFolder = '';
% 
sf_class   = {'description_sf_1y'};
%path_group = {[0 1 2] [0 1 2 3 4]};
path_group = {[0]};
bioNames    = {'ARR','PAC','PLV','PLI','H2','GC','sdDTF'};



for b = 1 : numel(bioNames)
    cfg.tbl2load{b} = fullfile(rootSummaryFolder,strcat('summary_tbl_',bioNames{b},'.mat'));
end

cfg.maxNout        = 0;
cfg.bioNames       = bioNames;
                


cfg.normTempRegExp = '1a_AED_stop\w*';

%cfg.seizOut2try    = {'1a_AED_stop\w*','1(a|b)\w*','\w*','(1(c|d)|2|3|4)\w*'};
cfg.seizOut2try    = {'1a_AED_stop\w*','1(a|b)\w*'};


cfg.path_label     = {'all'};

cfg.subj2rem       = {                                                        ...                                              
                      'RESP0448','RESP0480','RESP0482','RESP0484','RESP0497', ... % hfo trial
                      'RESP0500','RESP0519','RESP0537','RESP0542','RESP0556', ... % hfo trial
                      'RESP0566','RESP0572','RESP0604','RESP0631','RESP0636', ... % hfo trial
                      'RESP0640','RESP0644', ...                                  % hfo trial 
                      'RESP0353','RESP0623'  ...                                  % no last post recording
                     };

typeEPI = {'\w*','T','E'};  
% typeEPI = {'\w*'}; 

 out = [];    

 for e = 1: numel(typeEPI) 
     
     cfg.typeEPI = typeEPI{e};  
     
         for j = 1 : numel(sf_class)
                   
            cfg.sf_var      = sf_class{j}; 
            cfg.path_groups = path_group{j};

            %out{e,j} = risk_of_surgery(cfg);
            out{e,j}  = get_max_tbls_for_biomarker(cfg);
         end
    
 end
