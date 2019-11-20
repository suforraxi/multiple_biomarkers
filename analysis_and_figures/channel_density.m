function channel_density(cfg)
%close all

densityXChType     = [];
g_cardinality      = [];

f1 = figure;
f2 = figure;
subplot(4,2,1:2)
title(sprintf('Kolmogorov-Smirnov test Pre-resection \n resected channel distribution > Post-resection channel distribution'),'FontSize',16)

for t = 1 : numel(cfg.tblLoaded)
    subj_tbl = cfg.tblLoaded{t};
    [c_DxCT,vals,group_labels,g_cardinality,th_healthy] = density_of_channels_compared_to_normal(subj_tbl,cfg.prc_healthy); 
    densityXChType                           = [densityXChType; c_DxCT];
    figure(f1)
    %figure
    subplot(4,2,t)
    [h,p] = plot_distribution(cfg,vals,group_labels,cfg.bioName{t});
    switch cfg.bioName{t}
       
        case 'GC'     
            set(gca,'XLim',[0 0.005] );
        case 'sdDTF'  
            set(gca,'XLim',[0 0.0002] );
    end
    
    if(cfg.linePercentile)
        line([th_healthy th_healthy],[0 0.05],'Color','green','LineStyle','--')
        title(sprintf('%s th95: %.5f h:%i p:%.4f',cfg.bioName{t},th_healthy,h,p),'FontSize',13)
    else
        %title(sprintf('Kolmogorov-Smirnov test Pre-Resected > Post  p:%.4f',p),'FontSize',13)
    end
      figure(f2)
      subplot(2,4,t)
     
      preR_idx = strcmp(group_labels,'RES');
      post_idx = strcmp(group_labels,'Cured-PostRes');
      
      pre_post_vals  = [vals(preR_idx); vals(post_idx)];
    
      pre_post_label = [repmat({'PRE Resected'},1,sum(preR_idx)) repmat({'POST'},1,sum(post_idx))];
      viol = violinplot(pre_post_vals, pre_post_label);
      viol(1).ViolinColor = [0 0 1];
      viol(2).ViolinColor = [1 0 0];
      %title(out{te,1}.risk_str_T{1}.Row{idx2plot(i)}) 
      
      
      title(cfg.bioName{t},'FontSize',13)
     
      switch cfg.bioName{t}
          case {'ARR','PAC'}
              ylabel(sprintf('%s per channel',cfg.bioName{t}))
          otherwise
              ylabel(sprintf('strength per channel',cfg.bioName{t}))
      end
      
      xMarker = 1.5;
      MSize   = 14;
      ymax    = max(pre_post_vals);
     
      
      if(p <= cfg.alpha_level) % plot asteriks if significant
            plot(xMarker,ymax,'*','MarkerSize',MSize)
      end  
    
end



% comparison of the biomarker 'normal' vs 'potentially pathological' channels
% post-resection channels (not resected) of seizure free subjects are considered the 'normal' 
% (look the filter for definition of seizure free) 
% pre-resection channels are considered potentially pathological 
% pre-resection channels are divided in
% resected     - both channels in the bipolar derivation are in the resected area 
% cut          - one of the channels in the bipolar derivation is in the resected area
% not resected - both channels in the bipolar derivation are not in the resected area
% artefactual channels are not considered

% density_for_chType - percentage of channels with a value bigger than the
%                      threshold of the biomarker obtained from the 'normal' channels
% vals               - raw values organized for the different groups (RES/CUT/NRES/normal)
% group_labels       - label for the the vals

function [density_for_chType,vals,group_labels,g_cardinality,th_healthy] = density_of_channels_compared_to_normal(subj_tbl,prc_healthy) 

%find 'healthy tissue'
resected_idx  = strcmp(subj_tbl.resected,'RES');
cut_idx       = strcmp(subj_tbl.resected,'CUT');
nresected_idx = strcmp(subj_tbl.resected,'NRES');

art_free      = ones(size(subj_tbl,1),1);% subj_tbl.artefact == 0;
cured_def     = ~cellfun(@isempty,regexpi(subj_tbl.description_sf_1y,'1a_AED_stop'));
seiz_free_def = ~cellfun(@isempty,regexpi(subj_tbl.description_sf_1y,'1(a|b)\w*'));

pre_idx  = ~cellfun(@isempty,regexp(subj_tbl.sitName,'SITUATION1\w*'));
post_idx = ~pre_idx;


healthy_val = subj_tbl.biomarker(cured_def & post_idx & art_free);


th_healthy = prctile(healthy_val,prc_healthy);


epileptic_val_res    = subj_tbl.biomarker(seiz_free_def & pre_idx & art_free  & (resected_idx));
epileptic_val_cut    = subj_tbl.biomarker(seiz_free_def & pre_idx & art_free & cut_idx);
epileptic_val_nres   = subj_tbl.biomarker(seiz_free_def & pre_idx & art_free & nresected_idx);

vals   = [epileptic_val_cut; epileptic_val_res; epileptic_val_nres ; healthy_val];

group_labels = [repmat({'CUT'},1,length(epileptic_val_cut))  ...
    repmat({'RES'},1,length(epileptic_val_res)) ...
    repmat({'NRES'},1,length(epileptic_val_nres)) ...
    repmat({'Cured-PostRes'},1,length(healthy_val))];

ug_label           = unique(group_labels);
density_for_chType = zeros(1,numel(ug_label));
g_cardinality      = zeros(1,numel(ug_label));
for i = 1 : numel(ug_label)
    
    c_idx   = strcmp(group_labels,ug_label{i});
    c_group = vals(c_idx);
    
    density_for_chType(i) = sum((c_group > th_healthy))/length(c_group);
    g_cardinality(i)      = length(c_group);
end



% plot distribution of POST vs PRE_R

function [h,p] = plot_distribution(cfg,vals,group_labels,bioName)

nStep      = cfg.nStep;
nStepQ     = cfg.nStepQ;
normalize  = cfg.normalize;
showpoints = cfg.showpoints;

post_idx   = strcmp('Cured-PostRes',group_labels);
preR_idx   = strcmp('RES',group_labels);

post = vals(post_idx);
pre  = vals(preR_idx);

maxBio = max([pre ; post]);
minBio = min([pre ; post]);
step   = (maxBio - minBio)/nStep; 

th_post = max(post);


binCenters = 0 : step : maxBio+2*step;

[nC,binC] = hist(post,binCenters);
[nI,binI] = hist(pre,binCenters);

if(normalize)
    nC = nC/sum(nC);
    nI = nI/sum(nI);
end

if(cfg.spline)
    step       = (maxBio - minBio)/nStepQ;
    qpoints    = 0 : step :maxBio + 2*step;   
    cs_pre     = spline(binI ,nI ,qpoints);
    cs_post    = spline(binC ,nC ,qpoints);

else
    qpoints    = binCenters;   
    cs_pre     = nI;
    cs_post    = nC;
end

plot(qpoints,cs_post,'b')
hold
plot(qpoints,cs_pre,'r')

[h,p] = kstest2(post,pre,'Tail','larger');
if(showpoints)
    plot(binC,nC,'b*')
    plot(binI,nI,'ro')
end

if(normalize)
    ylabel('Channel density','FontSize',11)
else
    ylabel('# channels')
end
xlabel(bioName,'FontSize',14)


