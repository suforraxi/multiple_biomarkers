%% creating a summary table 

% name / age / gender / # electrodes / # resected electrodes / primary
% pathology / 1y seizure outcome

age_gender_F = '/home/matteo/Desktop/openclinica/TAB_age_gender_second_attempt_2019-10-28-095433974.tsv';
info_F       = '/home/matteo/Desktop/replicate_analysis/info/info.tsv';
bio_F        = '/home/matteo/Desktop/replicate_analysis/summary_tables/summary_tbl_ARR.txt';



age_gender_T = readtable(age_gender_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T       = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
subj_tbl     = readtable(bio_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);

primary_path = {'High Grade Tumor (WHO III + IV)',...
                    'Low Grade Tumor (WHO I + II)',   ...
                    'MTS',                            ...
                    'FCD',                            ...
                    'no abnormalities',               ...
                    'cavernoma',                      ...
                    'gliosis/scar',                   ...
                    'AVM',                            ...
                    'malformation cort. development', ...
                    'TuberoSclerosis'
                };
            
            
            
cfg.sf_var          = 'description_sf_1y';
cfg.sf_regexp       = '1(a|b)\w*';


                  
cfg.subj2rem       = {                                                        ...                                              
                          'RESP0448','RESP0480','RESP0482','RESP0484','RESP0497', ... % hfo trial
                          'RESP0500','RESP0519','RESP0537','RESP0542','RESP0556', ... % hfo trial
                          'RESP0566','RESP0572','RESP0604','RESP0631','RESP0636', ... % hfo trial
                          'RESP0640','RESP0644', ...                                  % hfo trial 
                          'RESP0353','RESP0623'  ...                                  % no last post recording
                          'RESP0218','RESP0301'  ...                                  % frequency sample not 2048
                          'RESP0118' ...                                              % no resected channels  
                          };

cfg.typeEPI         = '\w*';

cfg.pathology_group = 0; 


% select artefact free
subj_tbl = select_artfree(subj_tbl);


subj_tbl = select_seizoutcome(subj_tbl,cfg.sf_var,cfg.sf_regexp);

% remove subject by name
subj_tbl = rem_subj_from_tbl(cfg.subj2rem,subj_tbl);


% select pathology group
subj_tbl = select_primary_pathology_group(subj_tbl,cfg.pathology_group);

%filter for epilepsy type

subj_tbl = select_typeEPI(subj_tbl,cfg.typeEPI);


% total number of channels per subject

%idx2keep  = ~cellfun(@isempty,regexp(subj_tbl.sitName,'SITUATION1\w*'));
%post_idx = ~cellfun(@isempty,regexp(subj_tbl.sitName,'SITUATION[^1]\w*'));
subj_pre_tbl  = select_situations(subj_tbl,'SITUATION1\w*');
subj_post_tbl = select_situations(subj_tbl,'SITUATION[^1]\w');

[sG_pre,sIDs_pre] = findgroups(subj_pre_tbl.subjName);

isCut    = @(x) sum(strcmp(x,'CUT'));
isNR     = @(x) sum(strcmp(x,'NRES'));
isRes    = @(x) sum(strcmp(x,'RES'));

cutXsubj  = splitapply(isCut,subj_pre_tbl.resected,sG_pre);
resXsubj  = splitapply(isRes,subj_pre_tbl.resected,sG_pre);
nresXsubj = splitapply(isNR,subj_pre_tbl.resected,sG_pre);

preCH_tbl = cell2table(sIDs_pre,'VariableNames',{'subjID'});
preCH_tbl = [ preCH_tbl  ...
              array2table([cutXsubj nresXsubj resXsubj],'VariableNames',{'CUT','NResPre','Res'});
             ];

[sG_post,sIDs_post] = findgroups(subj_post_tbl.subjName);

nresPostXsubj = splitapply(isNR,subj_post_tbl.resected,sG_post);

postCH_tbl = cell2table(sIDs_post,'VariableNames',{'subjID'});
postCH_tbl = [ postCH_tbl  ...
              array2table(nresPostXsubj,'VariableNames',{'PostNR'});
             ];

% cured patients
cured_subj_tbl      = select_seizoutcome(subj_tbl,cfg.sf_var,'1A_AED_stop');
cured_subj_post_tbl = select_situations(cured_subj_tbl,'SITUATION[^1]\w');

[cured_sG_post,cured_sIDs_post] = findgroups(cured_subj_post_tbl.subjName);
cured_nresPostXsubj             = splitapply(isNR,cured_subj_post_tbl.resected,cured_sG_post);



info_T       = info_T(:,{'subjID','description_sf_1y','typeEPI','primary_path_class'});

age_gender_T = age_gender_T(:,{'StudySubjectID','Sex','DateOfBirth'});
age_gender_T.Properties.VariableNames = {'subjID','Gender','Age'};

current_year  = datetime(date).Year;
dateOFbirth   = age_gender_T.Age;
dateOFbirth   = dateOFbirth.Year;
newAge        = repmat(current_year,numel(dateOFbirth),1) - dateOFbirth;

age_gender_T.Age =  newAge;



characteristic_tbl = innerjoin(info_T,preCH_tbl,'Keys','subjID');
characteristic_tbl = innerjoin(characteristic_tbl,postCH_tbl,'Keys','subjID');


%varNames = {'subjID','sf_1y','typeEPI','CUT','NRes','Res'}

char_tbl = innerjoin(age_gender_T,characteristic_tbl,'Keys','subjID');

[~,idx_s] = sort(char_tbl.typeEPI);

char_tbl = char_tbl(idx_s,:);

char_tbl.description_sf_1y = lower(char_tbl.description_sf_1y);

% mean age
avg_age = mean(char_tbl.Age);
% number of male/female
nOfMale   = sum(strcmp('m',char_tbl.Gender));
nOfFemale = sum(strcmp('f',char_tbl.Gender));
% number of 1a_AED_stop / 1a_AED_eq / 1a_AED_low

[sG,typeSO] = findgroups(char_tbl.description_sf_1y);

countTypeSO = splitapply(@numel,char_tbl.description_sf_1y,sG);

% number of subjects per different pathology

[sG,typePath] = findgroups(char_tbl.primary_path_class);

countPriPath = splitapply(@numel,char_tbl.primary_path_class,sG);


% total number of pre-resected /pre noresected/pre cut / post channels
numChXtype = sum(char_tbl(:,{'CUT','NResPre','Res','PostNR'}).Variables);

% mean number of different type of channels across subject
avgChXtype = mean(char_tbl(:,{'CUT','NResPre','Res','PostNR'}).Variables);

sprintf('mean age: %.1f\nnumber of male/female: %i/%i\n',avg_age,nOfMale,nOfFemale)

summary_typeSO_tbl   = array2table(countTypeSO,'VariableNames',{'num'},'RowNames',typeSO);
summary_typePath_tbl = array2table(countPriPath,'VariableNames',{'num'},'RowNames',primary_path(typePath));
summary_chXtype_tbl  = array2table([numChXtype;avgChXtype],'VariableNames',{'CUT','NResPre','Res','PostNR'});

summary_typeSO_tbl
summary_typePath_tbl
summary_chXtype_tbl

prim_path_c     = char_tbl.primary_path_class;
prim_path_label = cell(size(prim_path_c));
for i = 1 : numel(prim_path_c)
    prim_path_label{i} = primary_path{prim_path_c(i)};
end

char_tbl.primary_path_class = prim_path_label;
