%
%  Compute results and figures separately for different biomarkers 
%  and combining them.
%
%  biomarkers (ARR / PAC / PLV / PLI / H2 / GC / sdDTF)
%
%  The analysis is divided in three main parts:
% 
% 1) Difference between distributions (pooling together all the available channels of different subjects) of 
%    pre-resection resected biomaker values in improved patients
%    and post-resection biomarker values in cured patients 
% 
% 2) Difference between distributions (maxim value of the biomarker per subjects) of 
%    pre-resection resected biomaker values in improved patients
%    and post-resection biomarker values in cured patients
% 
% 3) Combining all the biomarkers together 
%
%
% INPUT
% cfg is a struct with the following fields
% 
%       info_F                        - file name of the table information (see subj_info_F batch_compute_different_biomarkers)
%       bidsFolder                    - folder name where the raw data is
%                                       stored in BIDS format
%       rootInResFolder               - folder name where the results from
%                                       the biomarker computation are
%                                       stored. The stored struct of the result is
%
%                                                       outres.hdr          - containing the datasetName from which the
%                                                                             biomarker is computed
%                                                       outres.label        - cell array containing the channel names for which
%                                                                             the biomarker is computed 
%                                                       outres.time         - cell array containing the time for each trial in
%                                                                             which the original data is divided
%                                                       outres.fsample      - sample frequency of the data
%                                                       outres.sampleinfo   - sample information corresponding to the samples
%                                                                             in the original recordings from which each trial
%                                                                             computed
%                                                       outres.bio          - cell array with an entry per trial of the input
%                                                                             data. Each entry corresponds another cell array
%                                                                             with the computation of the biomarker per each
%                                                                             channel specified by label
%                                                                             (or biomaker statistic, like the strenght for bivariate/ multivariate methods)
%                                                       outres.extra        - struct with extra results depending on the
%                                                                             biomaker (see biomarker_wrapper specific function)
%                                                       outres.type         - name of the biomaker computed
%
%
%       rootSummaryFolder             -  folder name where to save the overall table for a biomarker   (see create_summary_table_main for the layout of the table)
%       poolingChannelFile            -  file name to save the group level analysis pooling the channel together
%       outFolderMaxComparison        -  folder name to save the figures of the maximum distributions comparison and the combining biomarker analysis        
%       subj2rem                      -  cell array of subject to remove
%                                        from the analysis i.e {'RESPXXX','RESPYYY'}
%       typeEPI                       -  cell array of regular expressions (see rem_subj_from_tbl )
%                                        indicating the type of epilepsy i.e. {'\w*'  'T'  'E'} (see select_typeEPI)
%       typeEPI_label                 -  cell array of labels for type of epilepsy i.e. {'Joint','Temporal',Extra-Temporal} 
%       sf_class                      -  cell array describing the seizure
%                                        outcome variable to use i.e {'description_sf_1y'} or {'description_sf_longest'} 
%       bioNames                      -  cell array of possible biomarker to compute i.e {'ARR'  'PAC'  'PLV'  'PLI'  'H2'  'GC'  'sdDTF'}
%       bioMarker2plotMaxDistribution -  array of indexes of the biomakers
%                                        for which to plot the violoin plot
%                                        of the comparison (i.e. [1 2 3 4 5 6 7]). Violin plot relative to the comparison of the biomarker maximum
%                                        distributions between pre-resection resected channels and post-resection channels (respectively in improved and cured patients)
%       alpha_level                   -  alpha level for the two sample one-sided Kolmogorov-Smirnov test
%       path_group                    -  cell array of integers (see
%                                        primary pathology class in create_summary_table_main) if 0 means all the pathologies 
%       path_group_label              -  cell array of labels corresposing to the primary pathology class considered
%       path_idx_of_interest          -  index to relative to the primary
%                                        pathology class to study (see get_max_tbls_for_biomarker)
%       seizOut2try                   -  cell array of regular expressions relative to which group of subjects to select i.e. cured or improved {'1a_AED_stop\w*'  '1(a|b)\w*'}
%       seizOut_idx                   -  indexes relative to seizOut2try 
%       seizOut_poolingCH             -  cell array with a regular
%                                        expression to select the subjects for whom run the group level
%                                        analysis pooling channels i.e. {'1(a|b)\w*'} for improved patients
%       removeNNchannels              -  1 to remove the channels labeled N-N in the name 0 to leave them 
%       
%       
function compute_results_and_save_figure(cfg)




% create tables from biomarkers struct (outres) 

% Result folder
rootInResFolder   = cfg.rootInResFolder;
rootSummaryFolder = cfg.rootSummaryFolder;

inFolder = dir(rootInResFolder);
inFolder = inFolder(~ismember({inFolder.name},{'.','..'}));


for d = 1 : numel(inFolder)

    folderName     = fullfile(inFolder(d).folder,inFolder(d).name);    
    cfg.resFolder  = folderName;

    cfg.fName2save = fullfile(rootSummaryFolder,strcat('summary_tbl_',inFolder(d).name,'.mat'));
    
    if(~isfile(cfg.fName2save)) % if the summary table is not present create it
        create_summary_table_main(cfg)
    end
end


sf_class   = cfg.sf_class;

path_group = cfg.path_group;

bioNames   = cfg.bioNames;


for b = 1 : numel(bioNames)
    cfg.tbl2load{b} = fullfile(rootSummaryFolder,strcat('summary_tbl_',bioNames{b},'.mat'));
end


cfg.bioNames       = cfg.bioNames;
                

typeEPI = cfg.typeEPI;  
  
% Difference between distributions (pooling together all the available channels of different subjects) of 
% pre-resection resected biomaker values in improved patients

pooling_channels(cfg)

% compute maximum table per type of epilepsy / seizure outcome group
% Each element of the cell out contains a struct with three dimensional cell biomaker X pathology group X seizure outcome group
% each cell contains a table with the following variables
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
out           = [];

 for e = 1: numel(typeEPI) 
     
     cfg.typeEPI = typeEPI{e};  
     
         for j = 1 : numel(sf_class)
                   
            cfg.sf_var      = sf_class{j}; 
            cfg.path_groups = path_group{j};

            out{e,j}        = get_max_tbls_for_biomarker(cfg);
         end
    
 end

 
% Difference between distributions (maximum value of the biomarker per subjects) of 
% pre-resection resected biomaker values in improved patients
% and post-resection biomarker values in cured patients 
 
 
cfg.path_idx    = cfg.path_idx_of_interest;
cfg.bioNames    = cfg.bioNames; 
cfg.idx2plot    = cfg.bioMarker2plotMaxDistribution; 
cfg.alpha_level = cfg.alpha_level;
cfg.outFolder   = cfg.outFolderMaxComparison;
cfg.Msize       = 14;
cfg.seizOutIdx  = [cfg.seizOut_idx];


% create output folder

if(~isfolder(cfg.outFolderMaxComparison))
    mkdir(cfg.outFolderMaxComparison);
end


compare_max_distributions(cfg,out);



% combining all biomarkers toghether as a 'cumulative' biomarker
combining_all_the_biomarkers(cfg,out)





% compare distributions pre-resection in improved patients vs post-resection cured pooling channels  
function pooling_channels(cfg)
figure

sf_class   = cfg.sf_class{1};
seizOut    = cfg.seizOut_poolingCH{1};
path_group = cfg.path_group{1};
typeEPI    = cfg.typeEPI{1};

for i = 1 :numel(cfg.tbl2load)
   
    load(cfg.tbl2load{i});
    
    
    % select artefact free
    subj_tbl = select_artfree(subj_tbl);


    subj_tbl = select_seizoutcome(subj_tbl,sf_class,seizOut);

    % remove subject by name
    subj_tbl = rem_subj_from_tbl(cfg.subj2rem,subj_tbl);

    
    % select pathology group
    subj_tbl = select_primary_pathology_group(subj_tbl,path_group);
    
    %filter for epilepsy type

    subj_tbl = select_typeEPI(subj_tbl,typeEPI);

    
    
    [vals,group_labels,~] = getValXclass_pre_post_R_NR(subj_tbl);


    preR_idx = strcmp(group_labels,'preR');
    post_idx = strcmp(group_labels,'postNR');
   
    [h,p]          = kstest2(vals(post_idx),vals(preR_idx),'Tail','larger');
    pre_post_vals  = [vals(preR_idx); vals(post_idx)];

    pre_post_label = [repmat({'PRE Resected'},1,sum(preR_idx)) repmat({'POST'},1,sum(post_idx))];
    
    subplot(2,4,i)
    
    viol = violinplot(pre_post_vals, pre_post_label);
    viol(1).ViolinColor = [0 0 1];
    viol(2).ViolinColor = [1 0 0];


    title(cfg.bioNames{i},'FontSize',13)

    switch cfg.bioNames{i}
      case {'ARR','PAC'}
          ylabel(sprintf('%s per channel',cfg.bioNames{i}))
      otherwise
          ylabel(sprintf('strength per channel',cfg.bioNames{i}))
    end

    xMarker = 1.5;
    MSize   = 14;
    ymax    = max(pre_post_vals);

    
    if(p <= cfg.alpha_level) % plot asteriks if significant
        plot(xMarker,ymax,'*','MarkerSize',MSize)
    end  


end


f = gcf;
f.PaperOrientation = 'landscape';
set(f, 'Position', get(0, 'Screensize'));
% create output folder
[outFolder,~,~] = fileparts(cfg.poolingChannelFile); 

if(~isfolder(outFolder))
    mkdir(outFolder);
end
print(cfg.poolingChannelFile ,'-dpng')
close all







% compare max distribution per subject with Kolmogorov-Smirnov test
function compare_max_distributions(cfg,out)

outFolder    = cfg.outFolder;
MSize        = cfg.Msize;
xMarker      = 1.5;
bioNames     = cfg.bioNames;
path_idx     = cfg.path_idx;
idx2plot     = cfg.idx2plot;
alpha_level  = cfg.alpha_level;
seizOutIdx   = cfg.seizOutIdx;

h = zeros(size(out,1),size(out{1}.max_T,1));
p = h;

for te = 1 : size(out,1)
   
    f2 = figure;
       
    for i = 1 : numel(idx2plot)
      
        post_val = out{te}.max_T{idx2plot(i),path_idx ,seizOutIdx(1)}.postNR_val;
        pre_val  = out{te}.max_T{idx2plot(i),path_idx ,seizOutIdx(end)}.preR_val;

        [h(te,i),p(te,i)] = kstest2(post_val,pre_val,'Tail','larger');

      
       subplot(1,numel(idx2plot),i)
       
       catNames = [repmat({'PRE resected'},size(pre_val),1); repmat({'POST'},size(post_val),1)] ;
       viol = violinplot([pre_val;post_val], catNames);
       viol(1).ViolinColor = [0 0 1];
       viol(2).ViolinColor = [1 0 0];
       
       c_bioName = bioNames{idx2plot(i)}; 
       title(c_bioName)
       
       ymax = max([pre_val;post_val]);
       
       
       switch c_bioName
          case {'ARR','PAC'}
              ylabel(sprintf('max %s per subject',c_bioName))
           
          otherwise
              ylabel(sprintf('max strength per subject',c_bioName))
       end
       
       if(round(p(te,i),2) <= alpha_level) % plot asteriks if significant
            plot(xMarker,ymax,'*','MarkerSize',MSize)
       end
    end   
    
    f2.PaperOrientation = 'landscape';
    set(f2, 'Position', get(0, 'Screensize'));
    print(sprintf('%smaxPreR_vs_Post_violinplot_%s',outFolder,cfg.typeEPI_label{te}),'-dpng')
    close(f2)
end





% using all biomarker together

% out cell with the different groups  
%     type of epilepsy            X seizure outcome definition 
% (Joint/Temporal/ExtraTemporal)  X (@1y / @longest)
% max_T cell with tables 
%      biomarker index          X primary pathology group              X seizure freedom class
% (ARR/PAC/PLV/PLI/H2/GC/sdDTF) X (all/High-Low Glioma/MST/FCD/ etc etc )  X (cured/improved/all)
% see ( get_max_tbls_for_biomarker )
function combining_all_the_biomarkers(cfg,out)
close all

% info table

info_F = cfg.info_F;
info_T = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T = info_T(:,{'subjID','primary_path_class','description_sf_1y','typeEPI'});


outFolder = cfg.outFolderMaxComparison;
fileNames = {'matrix_hits','subjhits_bar','atleast_bar'} ;
goi       = {'Joint (T+E)','Temporal','Extra-Temporal'};
bioNames  = cfg.bioNames;


sOutDef     = 1; % seizure outcome definition @1y 
toTest      = 2; % index of the subject table to test (improved subjects) 
pathIdx     = 1; % primary pathology index to use (all subjects)
idxCured    = 1; % index cured subjects

typeE       = [1 2 3]; % type of epilepsy index 1 joint (Temporal + Extra-Temporal)  / 2 Temporal / 3 Extra-Temporal 
noHitAnyHit = zeros(length(typeE),2);
f = [];
for te = typeE
    
    numBio   = size(out{te,sOutDef}.max_T,1);
   
    normTemp = zeros(numBio,1);
    % consider only the subjects with a value in as maximum, sometimes it
    % is not possible to compute the measure because there  were not enough
    % epochs (i.e. every epoch with at least one channel with artefacts)
    idxSujb2use   = ~isnan(out{te,sOutDef}.max_T{1,pathIdx,toTest}.preR_val);
    hits          = zeros(sum(idxSujb2use),numBio);


    for i = 1 : numBio

        % look for the maximum across all channels / all situations
        % post-resecion / all subjects cured after surgery (i.e. all 1A Engel who stopped medication after operation: 1A_AED_stop)
        normTemp(i,1) = max(out{te,sOutDef}.max_T{i,pathIdx,idxCured}.postNR_val);

        preRvals      = out{te,sOutDef}.max_T{i,pathIdx,toTest}.preR_val(idxSujb2use);
        % how many subjects in the pre-resection recordings have at least
        % one resected channel bigger than the threshold computed on post-resection
        % cured subjects 
        idx_aboveTh   = preRvals > normTemp(i);
        
        hits (:,i)     = idx_aboveTh;
    end

    subjNames   = out{te,sOutDef}.max_T{1,pathIdx,toTest}.subjName;
    idxSujb2use = ~isnan(out{te,sOutDef}.max_T{1,pathIdx,toTest}.preR_val);
    subjNames   = subjNames(idxSujb2use);


    bioHitXsubj = sum(hits,2); 
   


    subjName_T = cell2table(subjNames,'VariableNames',{'subjID'});
    hits_T     = array2table(hits,'VariableNames',bioNames);
    tot_hits_T = array2table(bioHitXsubj,'VariableNames',{'totBioHits'});

    cum_hits_T = [subjName_T hits_T tot_hits_T];

    usedSubj_T = innerjoin(cum_hits_T,info_T,'Keys','subjID');

    [~,idx_sorted] = sort(usedSubj_T.typeEPI);

    ordered_hits_T     = usedSubj_T(idx_sorted,bioNames);
    ordered_tot_hits_T = usedSubj_T(idx_sorted,{'totBioHits'});
    ordered_subjName_T = usedSubj_T(idx_sorted,{'subjID'});
    ordered_path_T     = usedSubj_T(idx_sorted,{'primary_path_class'});
    ordered_te_T       = usedSubj_T(idx_sorted,{'typeEPI'});


    hits_m       = ordered_hits_T.Variables;
    tot_hits_v   = ordered_tot_hits_T.Variables;
    path_v       = ordered_path_T.Variables;
    ord_subjName = ordered_subjName_T.subjID;

    m            = hits_m.*repmat(path_v,1,size(hits_m,2));


    for i = 1 : numel(ord_subjName)
        ord_subjName{i} = strcat(ord_subjName{i},'_',ordered_te_T.typeEPI{i});  
    end
    % matrix with primary pathology
    f(1) = figure;
    imagesc(m');

    primary_path = {'a) High Grade Tumor (WHO III + IV)',...
                    'b) Low Grade Tumor (WHO I + II)',   ...
                    'c) MTS',                            ...
                    'd) FCD',                            ...
                    'e) no abnormalities',               ...
                    'f) cavernoma',                      ...
                    'g) gliosis/scar',                   ...
                    'h) AVM',                            ...
                    'i) malformation cort. development', ...
                    'l) TuberoSclerosis'
                    };

    colorbar('Ticks',1:numel(primary_path),'TickLabels',primary_path)
    colormap('jet')
    xticks(1:numel(ord_subjName));
    xticklabels(ord_subjName)
    xtickangle(45)
    yticklabels(bioNames)
    title(goi{te},'FontSize',14)
    % how many biomarkers above the treshold X subject
    f(2) = figure;
    bar(tot_hits_v)
    xticks(1:numel(ord_subjName));
    xticklabels(ord_subjName)
    xtickangle(45)
    ylabel('Number of biomakers above threshold')
    xlabel('subjects')
    title(goi{te},'FontSize',14)

    cumBioHit = zeros(1,numBio);


    for i = 1 : numBio

        cumBioHit(i) = sum(bioHitXsubj >= i); 
    end

    noHits = sum(bioHitXsubj == 0);

    noHitAnyHit(te,1)   = noHits; 
    noHitAnyHit(te,end) = cumBioHit(1);
    
    cumBioHit = [noHits cumBioHit];
    numHits   = [0 1:numBio];
    % how many biomarkers above the treshold in terms of subjects
    f(3) = figure;
    numBioLabel = categorical({'0','at least 1',...
        'at least 2','at least 3',...
        'at least 4','at least 5',...
        'at least 6','at least 7'}); 
    bar(numBioLabel,cumBioHit)
    ylabel('Number of subjects')
    xlabel('Number of biomaker above threshold')
    title(goi{te},'FontSize',14)
    
    for i = 1 : numel(f)
        set(f(i), 'Position', get(0, 'Screensize'));
        saveas(f(i),fullfile(outFolder,strcat(fileNames{i},'_',goi{te})),'png')
    end
    
end
close all

figure;
typeE_label = {'Joint (T+E)','Temporal','Extra-Temporal'};
numBioLabel = categorical({'0','at least 1'}); 
for i = 1 : size(noHitAnyHit,1)
    
    subplot(1,size(noHitAnyHit,1),i)
    bar(numBioLabel,noHitAnyHit(i,:))
    ylabel('Number of subjects')
    xlabel('Number of biomakers above threshold')
    title(typeE_label{i},'FontSize',14)
end
f = gcf;
set(f, 'Position', get(0, 'Screensize'));
saveas(f,fullfile(outFolder,'noHitAnyHit_j_t_e'),'png')
close all

