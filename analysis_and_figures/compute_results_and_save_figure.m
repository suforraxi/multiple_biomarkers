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
%combining_all_the_biomarkers(cfg,out)

using_all_biomarkers_as_a_whole(cfg,out)


res_analysis.out                 = out;
res_analysis.seizOUT             = cfg.seizOut2try;
res_analysis.path_group_label    = cfg.path_group_label;
res_analysis.typeEPI             = cfg.typeEPI_label;
res_analysis.seizOutVariable     = cfg.sf_class;

save(fullfile(cfg.outFolderMaxComparison,'max_tables'),'res_analysis') % save max tables


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
       
       catNames = [repmat({'PRE resected'},size(pre_val),1 ); repmat({'POST'},size(post_val),1)] ;
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


% Results using combining all biomarkers together.
% Save the following figures
% - Three histograms corresponding to: the joint group (Temporal + Extra-Temporal), the Temporal group and Extra-Temporal group
%       * Histogram y-axis  number of biomarkers above the threshold (max biomarker across post-resection situations in cured subjects) x-axis subject name    
% - Histogram combining all biomarkers together: if at least one biomarker in the pre-resection resected channels is above the threshold it counts as a successful detection 
%   
%
function using_all_biomarkers_as_a_whole(cfg,out)

outFolder            = cfg.outFolderMaxComparison;
bioNames             = cfg.bioNames;
idx_improved         = 2; 
idx_cured            = 1;
idx_seizout_variable = 1;
idx_path_group       = 1;
idx_typeEPI          = [1 2 3]; % joint / temporal / extra-temporal

hitXtypeEPI = cell(1,length(idx_typeEPI));
numBio      = length(bioNames);
numSubjXbio = zeros(length(idx_typeEPI),numBio);
for te = idx_typeEPI  

    % check the maximum number of subject for whom there is a maximum values in
    % the pre-resection resected variable
    for i = 1 : numBio 
        numSubjXbio(te,i) = sum(~isnan(out{te,idx_seizout_variable}.max_T{i,idx_path_group,idx_improved}.preR_val));
    end

    [maxNumSubj, idx_Max_SubjNumber] = max(numSubjXbio(te,:));

    % number of subjects with at least one channel in the pre-resection that was resected and had a value bigger than the threshold computed on the cured  
    hits           = zeros(maxNumSubj,numBio);
    idx_subjNotNaN = ~isnan(out{te,idx_seizout_variable}.max_T{idx_Max_SubjNumber,idx_path_group,idx_improved}.preR_val);
    subjNameNotNaN = out{te,idx_seizout_variable}.max_T{idx_Max_SubjNumber,idx_path_group,idx_improved}.subjName(idx_subjNotNaN);
    hit_T          = array2table(hits,'VariableNames',bioNames,'RowNames',subjNameNotNaN);
    
    for i = 1 : numBio

        % look for the maximum across all channels / all situations
        % post-resecion / all subjects cured after surgery (i.e. all 1A Engel who stopped medication after operation: 1A_AED_stop)
        normThreshold = nanmax(out{te,idx_seizout_variable}.max_T{i,idx_path_group,idx_cured}.postNR_val);
        
        idxNotNaN     = ~isnan(out{te,idx_seizout_variable}.max_T{i,idx_path_group,idx_improved}.preR_val);
        
        subjNotNaN    = out{te,idx_seizout_variable}.max_T{i,idx_path_group,idx_improved}.subjName(idxNotNaN);
        preRvals      = out{te,idx_seizout_variable}.max_T{i,idx_path_group,idx_improved}.preR_val(idxNotNaN);
        % how many subjects in the pre-resection recordings have at least
        % one resected channel bigger than the threshold computed on post-resection
        % cured subjects 
        idx_aboveTh   = preRvals > normThreshold;
        subjBiggerTH  = subjNotNaN(idx_aboveTh);
        
        hit_T{subjBiggerTH,i} = 1;

       
    end
    
    hitXtypeEPI{te} = hit_T;
    
end



% order the joint set according to type of epilepsy

jointGroup_names = hitXtypeEPI{1}.Row;
idx_temporal     = 2;
for i = 1 : numel(jointGroup_names)
        if(any(strcmp(jointGroup_names{i},hitXtypeEPI{idx_temporal}.Row)))
            
            hitXtypeEPI{1}.Row{i} = strcat('T','_',hitXtypeEPI{1}.Row{i});
        else
            hitXtypeEPI{1}.Row{i} = strcat('E','_',hitXtypeEPI{1}.Row{i});
        end
end

[~, idx_sorted ]= sort(hitXtypeEPI{1}.Row);

hitXtypeEPI{1} = hitXtypeEPI{1}(idx_sorted,:); 

hitXtypeEPI;

% save per hits per subject 
goi = {'Joint (T+E)','Temporal','Extra-Temporal'};

for i = 1 : numel(hitXtypeEPI)

    numSubj     = size(hitXtypeEPI{i},1);
    totHitXsubj = sum(hitXtypeEPI{i}.Variables,2);
    f           = figure;
    bar(totHitXsubj)
    xticks(1:numSubj);
    xticklabels(hitXtypeEPI{i}.Row)
    a = gca;
    a.TickLabelInterpreter = 'none';
    xtickangle(45)
    ylabel('Number of biomakers above threshold')
    xlabel('subjects')
    title(goi{i},'FontSize',14);
    
    set(f, 'Position', get(0, 'Screensize'));
    saveas(f,fullfile(outFolder,strcat('hits_per_subjects','_',goi{i})),'png')
    close(f);
end

% combining biomarkers save cumulative hits (no hits / number of hits across subjects)
goi         = {'Joint (T+E)','Temporal','Extra-Temporal'};
numBioLabel = categorical({'0','at least 1'}); 
f   = figure;
for i = 1 : numel(hitXtypeEPI)

    totHitXsubj = sum(hitXtypeEPI{i}.Variables,2);
    
    subplot(1,numel(hitXtypeEPI),i)
    noHit            = sum(totHitXsubj == 0); 
    numberAtLeastOne = sum(totHitXsubj ~= 0);
    
    bar(numBioLabel,[noHit numberAtLeastOne]);
    ylabel('Number of subjects')
    xlabel('Number of biomakers above threshold')
    title(goi{i},'FontSize',14)
        
end
set(f, 'Position', get(0, 'Screensize'));
saveas(f,fullfile(outFolder,'combiningBio'),'png')
close(f);



%%
sum_T = [];
for i = 1 : numel(hitXtypeEPI)
    
    sum_T = [sum_T; varfun(@sum,hitXtypeEPI{i})];
    
end
sum_T.Properties.RowNames = {'Joint','Temporal','Extra-temporal'};

sum_T.Properties.VariableNames = replace(sum_T.Properties.VariableNames,'sum_','');

greaterThXBio = cell(size(numSubjXbio,1),size(numSubjXbio,2)); 
for te = 1 : size(numSubjXbio,1)
    for bio = 1 : size(numSubjXbio,2)
    
        greaterThXBio{te,bio} = strcat(num2str(sum_T{te,bio}),' / ',num2str(numSubjXbio(te,bio)));
    end
end

greaterThXBio_T = cell2table(greaterThXBio,'VariableNames',sum_T.Properties.VariableNames,'RowNames',sum_T.Properties.RowNames);

writetable(greaterThXBio_T,fullfile(outFolder,'typeEpiXBiomarker_table'),'FileType','text','Delimiter','tab','WriteRowNames',1)

