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

% INPUT 
% configuration structure cfg
%
% cfg.tbl2load     - cell with fileNames with the tables to load
% cfg.path_groups  - pathology group index to compute (see primary pathology class in create_summary_table_main)
% cfg.seizOut2try  - cell of regular expressions to define the seizure outcome groups where to test the global
%                    threshold
% cfg.typeEPI      - cell of regular expressions defining the epilepsy type (/w* = all , T = temporal, E = extra-temporal ) 
% cfg.subj2rem     - cell with subject codes to exclude from the analysis
%
% OUTPUT
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

%
% INPUT
% struct with the following fields
% cfg.ctbl2load        cell with file name of the table to load
% cfg.sf_var           variable to consider for the seizure outcome (description_sf_1y or description_sf_longest see create_summary_table_main)
% cfg.maxXcondRegExp   cell with a regular expression to define the seizure outcome groups where to test the global threshold
% cfg.subj2rem         cell with subject codes to exclude from the analysis
% cfg.pathology_group  integer relative to the pathology group (1:10 see create_summary_table_main) 
% cfg.typeEPI          cell of regular expression defining the epilepsy type (/w* = all , T = temporal, E = extra-temporal see create_summary_table_main) 
% 
% OUTPUT
% maxXcond_T        table with the following fields 
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

