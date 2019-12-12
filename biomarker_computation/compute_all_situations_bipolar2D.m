   
%
% Compute a biomarker for every subjects in BIDS folder
% Possible biomarkers:
%   AutoRegressive model Residual (ARR) (Geertsma 2017)
%   Phase Amplitude Coupling (PAC)
%   Phase Lag Index (PLI) (Stam 2007)
%   Phase Locking Value (PLV) (Mormann 2000)
%   H2 non linear correlation coefficient (Kalitzin 2006)
%   time-based Granger Causality (GC) (Lionel Barnett and Anil K. Seth, 2014 MVGC Toolbox)
%   Short-time direct Directed Transfer Function (sdDTF)(Mullen 2014 SIFT toolbox)
% INPUT
% cfgBatch.inDir_data    folder with raw data in BIDS format        
%
% cfgBatch.subj_info_F   file name of the table with information about the subjects (see batch_compute_different_biomarkers for the description of the table)   
%
% cfgBatch.outdir_combi  output folder where to save the biomarker
%                        computation for grid and strip combined (see batch_compute_different_biomarkers for the description of the output struct)
%
% cfgBatch.errorFile     error log file to write the failures during the
%                        computation

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
function compute_all_situations_bipolar2D(cfgBatch)


% directory to settings

inDir_data   = cfgBatch.inDir_data   ; 
subj_info_F  = cfgBatch.subj_info_F  ;

outdir_combi = cfgBatch.outdir_combi ;
errorFile    = cfgBatch.errorFile    ;


% obtain data (root directory of BIDS project)

subjList       = dir(fullfile(inDir_data,'sub*'));
subjNames_DB   = cell(numel(subjList),1);
for s = 1 : numel(subjList)  
    subjName        = regexp(subjList(s).name,'sub-(\w*)','tokens');
    subjNames_DB{s} = char(subjName{1});
end

% selection of  subjects 

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

% look if the format is specified in every recordings
for s = 1 : numel(subjList)
    
    sitDir = dir(fullfile(subjList(s).folder,subjList(s).name,'*SITUATION*'));
    
    if(~isempty(sitDir))

        for f = 1 :  numel(sitDir)
            sitFolder   = fullfile(sitDir(f).folder,sitDir(f).name,'ieeg');

            jfile2load  = dir(fullfile(sitFolder,'*ieeg.json'));

            jsonSideCar = loadjson(fullfile(jfile2load(1).folder,jfile2load(1).name));

            if (contains(jsonSideCar.iEEGElectrodeGroups,'Format;'))
                subjTODO{k,1}   = sitFolder;
                k = k + 1;
            end
        end
    end
end

% if already computed skip the data
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
     
   
    % cut the last epoch
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
    if(status == 0) % if the data was correctly imported
        
        ntrials = size(data.trial,2);
        % select last trial
        cfgLastEp        = [];
        cfgLastEp.trials = ntrials;

        data = ft_selectdata(cfgLastEp,data); 

        for f_el = 1 : numel(field_El)

           cfg = setfield(cfg,field_El{f_el},getfield(cfgBatch,field_El{f_el}));

        end
        
            % apply montage and compute biomarker
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


% Apply montage and compute biomarker
% cfg struct with the following fields
%
% datasetName  : file name (.vhdr) that will be analysed (see BIDS format)
% channelFile  : channel file name relative to the channel file in BIDS (see BIDS format) 
% annotFile    : annotation file we used a custom annotation file for ioECoG (see BIDS + https://github.com/suforraxi/ieeg_respect_bids )
% noArtefact   : field to remove it is not used anymore
% inDir_data   : input root directory for the data in BIDS
% subj_info_F  : file name containing a subject information table with variables specified in (batch_compute_different_biomarkers.m)
% cutLast      : 1 in order to cut the length amount of data (last minute in our case reduce propofol effect) of the recordings 
%                0 otherwise    
% trials       : fieldtrip field used to implement the cut (see ft_preprocessing) 
% length       : fieldtrip field used to implement the cut (see ft_preprocessing) 
% overlap      : fieldtrip field used to implement the cut (see ft_preprocessing) 
% cutTrials    : 1 in order to redefine the trials length 
%                0 otherwise
% trials_ct    : fieldtrip field used to implement the cut (see ft_redefinetrials) 
% length_ct    : fieldtrip field used to implement the cut (see ft_redefinetrials) 
% overlap_ct   : fieldtrip field used to implement the cut (see ft_redefinetrials) 
% deT_deM      : 1 to apply detred and demean (see ft_preprocesing)
%                0 otherwise               
% notch        : 1 to apply notch filter
%                0 otherwise   
% notchBS      : frequency interval where to apply the notch filter 
%                [low high] (see ft_preprocessing)
% montage      : 'bipolar_two_directions' apply the bipolar montage for
%                 grid 5x4 and strip 1x6 or 1x8 (custom function to compute montage see /montage/ folder)
% outdir_combi : folder name where to save the results (outres see below)
% errorFile    : file name where to save the failures
% epiBio       : biomaker name to compute it could be (ARR / PAC / PLV / PLI / H2 / GC / sdDTF )
%             
% extra field required, depending on the biomarker
%  ARR  (see wrapper ARR.m)
%           windowL : length in samples of the slinding window used to compute the ARR   
%  PAC  (see wrapper_PAC.m)    
%           lb: low frequency band boundaries [x y] used to estimate the phase   
%           hb: high frequency band boundaries [x y] used to estimate the
%               amplitude envelope
%  PLI  (see wrapper_PLI.m)
%           boi : [x y] frequency band boundaries to filter the signals before to compute PLI
%  PLV  (see wrapper_PLV.m)
%           boi : [x y] frequency band boundaries to filter the signals before to compute PLV
%  H2   (see wrapper_H2.m) 
%           boi : [x y] frequency band boundaries to filter the signals before to compute H2
%           T   : specify an array of delays to compute H2 (see external/h2delay.m)
%           n   : number of bins see (external/h2delay.m)
%  GC   (see wrapper_GCtime.m)
%           momax : maximal model order to try
%  sdDTF (see wrapper_sdDTF.m)
%           morder            :  model order 
%           freqBand          :  array specifying frequencies of interest
%           WindowLengthSec   :  window size for trial
%           WindowStepSizeSec :  step between windows
%
%
% output arguments
% status 0 the computation was ok 
%        1 error during the computation
% msg    error message in case of error
% outres The generic result is saved in a matlab structure with the following fields
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
           

function [status,msg,outres] = compute_montage_and_bio(cfg,data)
            
        status  = 0;
        msg    = [];
        outres = [];
        
        if(strcmp(cfg.montage,'bipolar_two_directions'))
            [status,msg,outres]    = compute_biomarker(cfg,data);
        end
        
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

