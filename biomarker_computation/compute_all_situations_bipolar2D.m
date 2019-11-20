
% Compute a biomarker for every subjects in BIDS folder
% Possible biomarkers:
% AutoRegressive model Residual (ARR) (Geertsma 2015)
% Phase Amplitude Coupling (PAC)
% Phase Lag Index (PLI) (Stam 2007)
% Phase Locking Value (PLV) (Mormann 2000)
% H2 non linear correlation coefficient (Kalitzin 2006)
% time-based Granger Causality (GC) (Lionel Barnett and Anil K. Seth, 2014 MVGC Toolbox)
% Short-time direct Directed Transfer Function (sdDTF)(Mullen 2014 SIFT toolbox)

% cfgBatch.inDir_data    BIDS folder with the data       
% cfgBatch.subj_info_F   table with information about the subjects (available/pathology/seizure outcome etc etc)   

% cfgBatch.outdir_combi  output folder where to save the biomarker
%                        computation for grid and strip combined
% cfgBatch.errorFile     error log file to write the failures

function compute_all_situations_bipolar2D(cfgBatch)


%% directory to settings

inDir_data   = cfgBatch.inDir_data   ; 
subj_info_F  = cfgBatch.subj_info_F  ;

outdir_combi = cfgBatch.outdir_combi ;
errorFile    = cfgBatch.errorFile    ;


%% obtain data (root directory of BIDS project)

subjList       = dir(fullfile(inDir_data,'sub*'));
subjNames_DB   = cell(numel(subjList),1);
for s = 1 : numel(subjList)  
    subjName        = regexp(subjList(s).name,'sub-(\w*)','tokens');
    subjNames_DB{s} = char(subjName{1});
end

%% selection of  subjects 

tbl_info_T      = readtable(subj_info_F, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);
tbl_info_T.Row  = tbl_info_T.subjID; 

% look for avaialable subjects
subj_selection_T = tbl_info_T(tbl_info_T.available == 1,:);
subj_selected    = subj_selection_T.subjID;

subj2use = false(numel(subjNames_DB),1);

% filter out not selected subjects
for s = 1:numel(subj_selected)
    
    subj2use = subj2use | strcmp(subj_selected{s},subjNames_DB);
end

subjList(~subj2use) = [];
subjTODO            = [];
k                   = 1;

% filter according to a specific format
for s = 1 : numel(subjList)
    
    sitDir = dir(fullfile(subjList(s).folder,subjList(s).name,'*SITUATION*'));
    
    if(~isempty(sitDir))

        for f = 1 :  numel(sitDir)
            sitFolder   = fullfile(sitDir(f).folder,sitDir(f).name,'ieeg');

            jfile2load  = dir(fullfile(sitFolder,'*ieeg.json'));

            jsonSideCar = loadjson(fullfile(jfile2load(1).folder,jfile2load(1).name));

            if (contains(jsonSideCar.iEEGElectrodeGroups,'Format;'))%'Format;Gr[4x5]'
                subjTODO{k,1}   = sitFolder;
                k = k + 1;
            end
        end
    end
end

% filter already computed
subjTODO;
bioDir   =  cfgBatch.epiBio;
idx2Skip =  false(size(subjTODO));

for i = 1 : numel(subjTODO)
       
       fileName = dir(fullfile(subjTODO{i},'*.vhdr'));
       [~,file2check,~] = fileparts(fullfile(fileName(1).folder,fileName(1).name));
       
       fileRes = (fullfile(outdir_combi,bioDir,strcat(file2check,'.mat')));
       if(isfile(fileRes))
            idx2Skip(i) = 1;
       end
end

subjTODO(idx2Skip) = [];

% create all the output directories and error files
[pathError,~,~] = fileparts(errorFile);
if(~isfolder(pathError)) 
    mkdir(pathError);
end
  
fid_bip = fopen(errorFile,'a+');


if(~isfolder(fullfile(outdir_combi,bioDir)))
    mkdir(fullfile(outdir_combi,bioDir));
end

% compute selected biomarker for every subjects

field_El = fields(cfgBatch);

for s = 1 : numel(subjTODO)

    [status,msg,cfg,data] = importBidsData(subjTODO{s});
     
   
    %% cut the last epoch
    if(cfgBatch.cutLast == 1)
        try
            cfgReTrials.trials  = cfgBatch.trials;
            cfgReTrials.length  = cfgBatch.length; %seconds of new trials
            cfgReTrials.overlap = cfgBatch.overlap;

            data = ft_redefinetrial(cfgReTrials,data);
        catch ME
            status = 1;
            msg = sprintf('%s err:%s --func:%s',subjTODO{s},ME.message,ME.identifier);
            data = [];
        end
    end
    if(status == 0)
        ntrials = size(data.trial,2);
        % select last trial
        cfgLastEp        = [];
        cfgLastEp.trials = ntrials;

        data = ft_selectdata(cfgLastEp,data); 

        for f_el = 1 : numel(field_El)

           cfg = setfield(cfg,field_El{f_el},getfield(cfgBatch,field_El{f_el}));

        end
        

            [status,msg,outres] = compute_montage_and_bio(cfg,data);

            if(status== 0 && ~isempty(outres))
                [~,file2save,~] = fileparts(cfg.datasetName);    
                save(fullfile(outdir_combi,bioDir,(file2save)),'outres'); 
            else % status
                fprintf(fid_bip,'%s',msg);
            end 

    else % status
            fprintf(fid_bip,'%s',msg);
    end
    
end % subjTODO

fclose(fid_bip);


%% montage & biomarker
function [status,msg,outres] = compute_montage_and_bio(cfg,data)
            
        status  = 0;
        msg    = [];
        outres = [];
        
        if(strcmp(cfg.montage,'bipolar_two_directions'))
            [status,msg,outres]    = compute_biomarker(cfg,data);
        end
        
% generic biomarker computation for strip and grid together 
% Note: removing bad channels from the bipolar montage (N-N)
function [status,msg,outres] = compute_biomarker(cfg,data)
          
status = 0;
msg    = [];
try
    outres = [];
   
   
     m_data_grid  = compute_montage_grid(cfg,data);
     
     m_data_strip = compute_montage_strip(cfg,data);
     
     m_data  = merge_dataset( m_data_grid,m_data_strip);
     
     
    if(~isempty(m_data))
        outres = biomarker_wrapper(cfg,m_data);
    end
    
catch ME
   
   status = 1;
   msg = sprintf('%s err:%s --func:%s',cfg.datasetName,ME.message,ME.stack(1).name);
   outres = [];
end
function m_data = compute_montage_grid(cfg,data)
    m_data = [];
    % apply a montage
    switch cfg.montage 
        case 'bipolar_two_directions'
    
            mX_data = apply_montage2data(data,@create_bipolar_montage_5by4X);
            mY_data = apply_montage2data(data,@create_bipolar_montage_5by4Y);
            m_data  = merge_dataset(mX_data,mY_data);
    end 

function m_data = compute_montage_strip(cfg,data)
    m_data = [];
    % apply a montage
    switch cfg.montage 
        case 'bipolar_two_directions'
            m_data = apply_montage2data_strip(data,@create_bipolar_montage_strip);
    end 

