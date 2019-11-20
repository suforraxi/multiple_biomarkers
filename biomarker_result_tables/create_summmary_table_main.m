% Create a general table with all the biomarker values per channel per situation 
% per subject, with resected information, seizure outcome and format

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
% have the same sitName/ subjName / paht_desc 


function create_summmary_table_main(cfg)

subj_tbl = create_table(cfg);
chName   = subj_tbl.chName;

fName2save = cfg.fName2save;


% remove bad channel N-N or not existing channel N-N 
if(cfg.removeNNchannels)
    aux      = regexp(chName,'\w*N-N\w*','match');
    idx2keep = cellfun(@isempty,aux);
    subj_tbl = subj_tbl(idx2keep,:);
end

save(fName2save,'subj_tbl');
writetable(subj_tbl,replace(fName2save,'.mat','.txt'),'Delimiter','tab');



function subj_tbl = create_table(cfg)

bidsFolder = cfg.bidsFolder;
resFolder  = cfg.resFolder;
info_F     = cfg.info_F;

% subject general info (imported correctly etc etc)

info_T = readtable(info_F,'Delimiter','\t','FileType','text','ReadRowNames',1,'ReadVariableNames',1);

% Information about the computation of the biomaker per situation
% for all the subjects in bids format check if the biomarker is computed.
% Check if the each subject has the biomarker computed for all the
% available situations. 
% Check if at least each subject has the biomarker computeed
% for at least one situation in the pre-resection condition and at least
% one in the post-resection condition.
% The output of the function is a table with RESPECT number as Rows
% and two varibles 
% all_completed ( 1 if the biomarker was computed for all the available situations, 0 otherwise )
% anyPre_anyPost ( 1 if the biomarker was computed for at least one situation pre an post, 0 otherwise )

cfg2.bidsDir = bidsFolder;
cfg2.resDir  = resFolder;

anypre_anypost_T = collect_info_about_situations_main(cfg2);

% The requisite is that there is a least one pre and one post situation per
% subject

idx_subj2use = anypre_anypost_T.anyPre_anyPost;

% format information
cfg3.bidsDir  = bidsFolder;
cfg3.tleInfoF = info_F;
format_T      = format_info_main(cfg3);

subjList = anypre_anypost_T.Properties.RowNames(logical(idx_subj2use));

subj_tbl = [];
for i = 1 : numel(subjList)

    c_subjName = subjList{i};
    
    idx_sits    = strcmp(c_subjName,format_T.subjNameXsit);
    
    c_format_T  = format_T(idx_sits,:);
    
    % consider only situation for the biomarker was computed
    idx2keep = true(size(c_format_T,1),1);
    for j = 1 : size(c_format_T,1) 
        
        c_fName = strcat(c_format_T.nameXsit{j},'.mat');
        if(~isfile(fullfile(resFolder,c_fName)))
            idx2keep(j) = 0;
        end
    end
    
    c_format_T = c_format_T(idx2keep,:);
    
    aux = regexp(c_format_T.nameXsit,'ses-(?<sitname>\w*)_','names');
    
    c_sitname = {};
    for j = 1 : numel(aux)
        c_sitname{j,1} = aux{j}.sitname;
    end
    
    [pre,inter,post] = find_pre_int_post(c_sitname);
    
    
    vars = [];
    % take only pre and post
    c_format_T       = c_format_T(pre | post,:);
    
    
    nSits            = size(c_format_T,1);
    
    desc_sf_1y       = info_T{c_subjName,'description_sf_1y'};
    desc_sf_longest  = info_T{c_subjName,'description_sf_longest'};
    path_desc        = info_T{c_subjName,'primary_path_class'};
    typeEPI          = info_T{c_subjName,'typeEPI'};
    HFOstudy         = info_T{c_subjName,'HFOstudy'};
    
    sitDesXsubj = [];
    for j = 1 : nSits
        format_entry = c_format_T(j,:); 
        
        sitDes_T     = get_sit_descriptor(bidsFolder,resFolder,...
                                          format_entry,...
                                          desc_sf_1y,...
                                          desc_sf_longest,...
                                          path_desc,...
                                          typeEPI,...
                                          HFOstudy);
        
        sitDesXsubj  = [sitDesXsubj ; sitDes_T];             
    end
    
    subj_tbl = [subj_tbl; sitDesXsubj];

end



function sitDes_T = get_sit_descriptor(bidsFolder,resFolder,format_T,desc_sf_1y,desc_sf_longest,path_desc,typeEPI,HFOstudy)

sitName_F = format_T.nameXsit{1};

load(fullfile(resFolder,strcat(sitName_F,'.mat')));

fileName_v = sitName_F;

aux        = regexpi(sitName_F,'\w*-(RESP\w*)_\w*','tokens');
subjName_v = aux{1}{1};

aux        = regexpi(sitName_F,'\w*-(SITUA\w*)_\w*','tokens');
sitName_v  = aux{1}{1};

format_v   = format_T.formatXsit{1};
nGrid_v    = format_T.nGridXsit(1);
nStrip_v   = format_T.nStripXsit(1);

description_sf_1y      = desc_sf_1y;
description_sf_longest = desc_sf_longest;



[ res_channel, artefact_T] = get_metadata(bidsFolder,sitName_F);

nCh                         = numel(outres.label);
[idxArtefact,idx_art_trial] = find_marked_artefacts(outres,artefact_T);

switch(outres.type)
    case {'ARR','PAC','PLV','PLI','H2','sdDTF'} 
        biores      = nanmean(abs(cell2mat(outres.bio)),2);
    case {'GC'}
        biores      = cell2mat(outres.bio);
end

biom_v = zeros(nCh,1);

for i = 1 : nCh  
    
    chName_v{i,1}   = outres.label{i};
    resected_v{i,1} = get_resected_label(chName_v{i,1},res_channel);
    
    artefact_v(i,1) = idxArtefact(i);

    if(~isempty(biores) && length(biores) == nCh)
        
        biom_v(i,1)     = biores(i);
       
    else
        biom_v(i,1)     = NaN;
    end
end

varNames = {'fName','subjName','sitName','chName','format','nGrid',...
            'nStrip','resected','artefact','biomarker','description_sf_1y','description_sf_longest','path_desc',...
            'typeEPI','HFOstudy'};
sitDes_T = table(repmat({fileName_v},nCh,1),   ...
                 repmat({subjName_v},nCh,1),   ...   
                 repmat({sitName_v},nCh,1),    ...
                 chName_v,                     ...
                 repmat({format_v},nCh,1),     ...
                 repmat(nGrid_v,nCh,1),        ...
                 repmat(nStrip_v,nCh,1),       ...
                 resected_v,                   ...
                 artefact_v,                   ...
                 biom_v,                       ...
                 repmat(description_sf_1y,nCh,1),  ...
                 repmat(description_sf_longest,nCh,1),...
                 repmat(path_desc,nCh,1),      ...
                 repmat(typeEPI,nCh,1),        ...
                 repmat(HFOstudy,nCh,1),       ...
                 'VariableNames',varNames      ...
             );
sitDes_T;

function res_notR_cut = get_resected_label(chName,res_channel)



ch = regexpi(chName,'(\w*\d+)[N]?-[N]?(\w*\d+)','tokens');

ch1 = ch{1}{1}; 
ch2 = ch{1}{end};


res_notR_cut = 'NRES';
if( any(strcmp(ch1,res_channel)) || any(strcmp(ch2,res_channel)) )
    if(any(strcmp(ch1,res_channel)) && any(strcmp(ch2,res_channel)))
        res_notR_cut = 'RES';
    else
        res_notR_cut = 'CUT';
    end
end

function [idxArtefact ,idx_art_trial] = find_marked_artefacts(outres,artefact_T)


s = outres.sampleinfo;


artefact_T;
nArtefacts = numel(artefact_T.type);
ch2remove  = [];
k          = 1;

idx_art_trial = zeros(1,size(s,1));
for i = 1 : nArtefacts
    for e = 1 : size(s,1)

        b_art   = artefact_T.start(i);
        end_art = artefact_T.stop(i);
        art_ch  = artefact_T.channel{i};

        if( ~( (s(e,1) > b_art && s(1)> end_art) || ( s(e,end) < b_art && s(e,end) < end_art) ) )
            if ( s(e,1) > b_art && s(e,end) > end_art )
                ch2remove{k}  = art_ch;
                k             = k+1; 
                idx_art_trial(e) = 1;
            end
            if( s(e,1) < b_art && s(e,end) < end_art  )
                ch2remove{k}  = art_ch;
                k             = k+1;
                idx_art_trial(e) = 1;
            end
            if( s(e,1) < b_art && s(e,end) > end_art  )
                ch2remove{k}  = art_ch;
                k             = k+1;
                idx_art_trial(e) = 1;
            end
            if( s(e,1) > b_art && s(e,end) < end_art  )
                ch2remove{k}  = art_ch;
                k             = k+1;
                idx_art_trial(e) = 1;
            end
        end
    end
end
idx_art_trial = logical(idx_art_trial);
ch2remove    = unique(ch2remove);
idxArtefact  = zeros(numel(outres.label),1);
if(~isempty(ch2remove))
    c_pattern = [];
    for i = 1 : numel(ch2remove)

        if(i==1)
            c_pattern = ['\w*' ch2remove{i} '\w*'];
        else
            c_pattern = [ c_pattern '|' '\w*' ch2remove{i} '\w*'];
        end

    end
        c_pattern = ['(' c_pattern ')'];

       aux      = regexpi(outres.label,c_pattern);
       idxArtefact = ~cellfun(@isempty,aux);
       
        
end









