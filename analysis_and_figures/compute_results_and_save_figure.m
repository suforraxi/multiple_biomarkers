%% 
%  Figures for manuscript DOI:  
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
function compute_results_and_save_figure(cfg)




%% create tables from biomarkers struct (outres) 

% Result folder
rootInResFolder   = cfg.rootInResFolder;
rootSummaryFolder = cfg.rootSummaryFolder;

inFolder = dir(rootInResFolder);
inFolder = inFolder(~ismember({inFolder.name},{'.','..'}));


for d = 1 : numel(inFolder)

    folderName     = fullfile(inFolder(d).folder,inFolder(d).name);    
    cfg.resFolder  = folderName;

    cfg.fName2save = fullfile(rootSummaryFolder,strcat('summary_tbl_',inFolder(d).name,'.mat'));
    
    create_summary_table_main(cfg)

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
out = [];  
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

compare_max_distributions(cfg,out);


%% using all biomarker together

% out cell with the different groups  
%     type of epilepsy            X seizure outcome definition 
% (Joint/Temporal/ExtraTemporal)  X (@1y / @longest)
% max_T cell with tables 
%      biomarker index          X primary pathology group              X seizure freedom class
% (ARR/PAC/PLV/PLI/H2/GC/sdDTF) X (all/High-Low Glioma/MST/FCD/ etc etc )  X (cured/improved/all)

close all

% info table

info_F = cfg.info_F;
info_T = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T = info_T(:,{'subjID','primary_path_class','description_sf_1y','typeEPI'});


outFolder = cfg.outFolderMaxComparison;
fileNames = {'matrix_hits','subjhits_bar','atleast_bar'} ;
goi       = {'Joint (T+E)','Temporal','Extra-Temporal'};


sOutDef     = 1; % seizure outcome definition @1y 
toTest      = 2; % index of the subject table to test (improved subject) 
pathIdx     = 1; % primary pathology index to use (all subject)
idxCured    = 1; % index cured subjects

typeE       = [1 2 3]; % type of epilepsy index 1 joint / 2 Temporal / 3 Extra-Temporal 
noHitAnyHit = zeros(length(typeE),2);
f = [];
for te = typeE
    
    numBio   = size(out{te,sOutDef}.max_T,1);
    numSubj  = size(out{te,sOutDef}.max_T{1,pathIdx,toTest},1);
    normTemp = zeros(numBio,1);
    
    idxSujb2use   = ~isnan(out{te,sOutDef}.max_T{1,pathIdx,toTest}.preR_val);
    hits          = zeros(sum(idxSujb2use),numBio);% one subject does not have resected channels for that situation


    for i = 1 : numBio


        normTemp(i,1) = max(out{te,sOutDef}.max_T{i,pathIdx,idxCured}.postNR_val);

        preRvals      = out{te,sOutDef}.max_T{i,pathIdx,toTest}.preR_val(idxSujb2use);
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
    print(sprintf('%smaxPreR_vs_Post_violinplot_%i',outFolder,te),'-dpng')
    close(f2)
end