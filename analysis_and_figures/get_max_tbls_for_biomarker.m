% extract maximum table from subject table (subj_tbl)
% 
% subject table is a table with the following variables:
%
% fName                 : filename of the result file corresponding to the situation and
%                         biomarker computed
% subjName              : coded subject name (RESPect name of the subject)
% sitName               : situation name (corresponding to the arrangement of the grid on the brain)
% chName                : recording site for which the biomarker value is computed
% format                : type of grid and strip used (i.e. 5x4, 4x8, 1x8, 1x6 )
% nGrid                 : number of grid used for the situation (sitName)
% nStrip                : number of strip used for the situation (sitName)
%
% resected              : RES  if the chName was resected 
%                         NRES if the chName was not resected
%                         CUT  (in bipolar montage)if the one chName of the bipolar derivarion was resected and the other not 
%
% artefact              : 1 if the chName contains some artefact 0 no artefact
% biomarker             : value of the computed biomarker for the chName
% description_sf_1y     : seizure freedom outcome one year after surgery (Engel class + medication
%                         level after surgery) 
%                        i.e. 
%                          1A_AED_eq   if the subject is Engel 1A class with
%                                      same medication after surgery
%                          1A_AED_stop if the subject is Engel 1A class who
%                                      stopped medication after surgery
%                          1A_AED_low  if the subject is Engel 1A class who 
%                                      decreased medication after surgery
%                 
% description_sf_longest : longest seizure freedom outcome reported after surgery  
%                         (Engel class + medication level after surgery) 
%
% paht_desc              : pathology descriptor, integer representing a
%                          primary pathology class
%                          1  High Grade Tumor (WHO III + IV)
%                          2  Low Grade Tumor (WHO I + II)
%                          3  MTS
%                          4  FCD
%                          5  no abnormalities
%                          6  cavernoma
%                          7  gliosis/scar
%                          8  AVM
%                          9  malformation cortical development
%                          10 TuberoSclerosis
%
% typeEPI                : type of epilepsy for the subject
%                          T Temporal
%                          E Extra Temporal
%
% HFOstudy               : 1 if the subject was in the HFO trial 
%                          0 otherwise

% NB there is redundancy in the information of the table since the table is built to
% have the main key variable as chName. 
% This means that all channel names in a situation will
% have the same sitName/ subjName / paht_desc / typeEPI etc etc.
% 

% input configuration structure cfg
%
% cfg.tbl2load     - filename with the table to load
% cfg.path_groups  - pathology group index to compute (see primary patology class)
% cfg.seizOut2try  - cell of regular expression to define the seizure outcome groups where to test the global
%                    threshold
% cfg.typeEPI      - cell of regular expression defining the epilepsy type (/w* = all , T = temporal, E = extra-temporal ) 
% cfg.subj2rem     - cell with subject codes to exclude from the analysis

% output struct with the following field
% out.max_T       -  three dimensional cell with a table for each cell relative to    
%              
%      biomarker index          X  primary pathology group                 X seizure freedom class
% (ARR/PAC/PLV/PLI/H2/GC/sdDTF) X (all/High-Low Glioma/MST/FCD/ etc etc )  X (cured/improved/)
%           
%                   each table has the following variables 
%
%                   subjName     - coded subject name
%                   postNR_sit   - situation name from post-resection where
%                                  the maximum value of the biomarker is found 
%                   post_chName  - channel name where the the maximum value of the biomarker is found
%                   postNR_val   - biomarker value corresponding to the maximum post-resection 
%                   preR_sit     - situation name from pre-resection recordings where
%                                  the maximum value of the biomarker is found across resected channels 
%                   preR_chName  - channel name where the the maximum value
%                                  of the biomarker is found across resected channels 
%                   preR_val     - biomarker value corresponding to the
%                                  maximum across pre-resection resected channels
%                   preNR_sit    - situation name from pre-resection recordings where
%                                  the maximum value of the biomarker is found across not resected channels
%                   preNR_chName - channel name where the the maximum value
%                                  of the biomarker is found across not resected channels
%                   preNR_val    - biomarker value corresponding to the
%                                  maximum across pre-resection not resected channels

function out = get_max_tbls_for_biomarker(cfg)

max_T = [];

for i = 1 : numel(cfg.tbl2load)
    
    cfg.ctbl2load = cfg.tbl2load{i};
    
    for pg = 1 : numel(cfg.path_groups) 
        
        cfg.pathology_group = cfg.path_groups(pg);
        
        for j = 1 : numel(cfg.seizOut2try)

            cfg.maxXcondRegExp = cfg.seizOut2try{j};

            max_T{i,pg,j} = get_max_tbls(cfg);
            
        end
    end
end

out.max_T = max_T;

function [maxXcond_T] = get_max_tbls(cfg)

condName = {'postNR','preR','preNR'};

subj_max   = [];


varNames = {'subjName','postNR_sit','post_chName','postNR_val',...
            'preR_sit','preR_chName','preR_val',...
            'preNR_sit','preNR_chName','preNR_val'};


% load table with biomarker info
load(cfg.ctbl2load)    

% select artefact free
subj_tbl = select_artfree(subj_tbl);

% select seizure outcome group
subj_tbl = select_seizoutcome(subj_tbl,cfg.sf_var,cfg.maxXcondRegExp);

% remove subject by name
subj_tbl = rem_subj_from_tbl(cfg.subj2rem,subj_tbl);

% subject selection based on the primary pathology group
subj_tbl = select_primary_pathology_group(subj_tbl,cfg.pathology_group);

%filter for epilepsy type

subj_tbl = select_typeEPI(subj_tbl,cfg.typeEPI);



for c = 1 : numel(condName)  
    subj_max{c} = get_MaxConditionTbl_pre_post_R_NR(subj_tbl,condName{c});
end

maxXcond_T = innerjoin(subj_max{1},subj_max{2},'Key','subjName');
maxXcond_T = innerjoin(maxXcond_T,subj_max{3},'Key','subjName');

maxXcond_T.Properties.VariableNames = varNames;     
    
 
% get max table per condition express in the regular expression condRegExp
function [aux_tbl] = get_MaxConditionTbl_pre_post_R_NR(subjT,condRegExp)

[subjG,ids] = findgroups(subjT.subjName);

subj_max  = zeros(numel(ids),1);

aux_tbl = [];

for s = 1 : numel(ids)
    
    idx_cs = subjG == s;
    
    c_subj_tbl = subjT(idx_cs,:);
    
    [vals,G_idx,new_tbl] = getValXclass_pre_post_R_NR(c_subj_tbl);
    
    ii = ~cellfun(@isempty,regexpi(G_idx,condRegExp,'match'));
    
    numNan = sum(isnan(vals(ii)));
    
    if(sum(ii) == 0 || numNan == length(vals(ii)))
        
        subj_max(s) = NaN;
        c_tbl = new_tbl(1,:);
        c_tbl.biomarker = NaN;
        c_tbl.sitName   = {'none'};
        c_tbl.chName    = {'none'};
        aux_tbl  = [aux_tbl; c_tbl];
        
    else
    
       [subj_max(s), max_idx] = max(vals(ii));
    
       c_tbl    = new_tbl(ii,:);
       aux_tbl  = [aux_tbl; c_tbl(max_idx,:)];
    end
end

% Extract the biomarker values according to the class (RES/CUT/NRES pre & post)
% and rename them in preR / preNR / postNR 

function [vals,G_idx,new_tbl] = getValXclass_pre_post_R_NR(subj_tbl)

resected_idx  = strcmp(subj_tbl.resected,'RES');
cut_idx       = strcmp(subj_tbl.resected,'CUT');
nresected_idx = strcmp(subj_tbl.resected,'NRES');

pre_idx  = ~cellfun(@isempty,regexp(subj_tbl.sitName,'SITUATION1\w*'));
post_idx = ~cellfun(@isempty,regexp(subj_tbl.sitName,'SITUATION[^1]\w*'));


preR   = subj_tbl.biomarker(pre_idx  & resected_idx );
preNR  = subj_tbl.biomarker(pre_idx  & nresected_idx);
postNR = subj_tbl.biomarker(post_idx & nresected_idx);

chName_preR   = subj_tbl.chName(pre_idx  & resected_idx);
chName_preNR  = subj_tbl.chName(pre_idx  & nresected_idx);
chName_postNR = subj_tbl.chName(post_idx & nresected_idx);


vals = [preR ;preNR ; postNR];
G_idx = [repmat({'preR'},size(preR,1),1) ; repmat({'preNR'},size(preNR,1),1) ; repmat({'postNR'},size(postNR,1),1)];

sName_preR   = subj_tbl.sitName(pre_idx  & resected_idx );
sName_preNR  = subj_tbl.sitName(pre_idx  & nresected_idx);
sName_postNR = subj_tbl.sitName(post_idx & nresected_idx);

subName_preR   = subj_tbl.subjName(pre_idx & resected_idx);
subName_preNR  = subj_tbl.subjName(pre_idx & nresected_idx);
subName_postNR = subj_tbl.subjName(post_idx& nresected_idx);


val_tbl  = array2table([ preR ; preNR ; postNR],'VariableNames',{'biomarker'});
ch_tbl   = array2table([ chName_preR  ; chName_preNR ; chName_postNR],'VariableNames',{'chName'});
sit_tbl  = array2table([ sName_preR ; sName_preNR ; sName_postNR],'VariableNames',{'sitName'});
sub_tbl  = array2table([ subName_preR ; subName_preNR ; subName_postNR],'VariableNames',{'subjName'});



new_tbl  = [sub_tbl sit_tbl ch_tbl val_tbl];
