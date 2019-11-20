%% group stats
close all
path_idx = 1;
h = zeros(size(out,1),size(out{1,1}.max_risk_T,1));
p = h;
for te = 1 : size(out,1)
    f1 = figure;
    f2 = figure;
    
    idx2plot = [1 2 3 4 5 6 7 ];
    
    for i = 1 : numel(idx2plot) %1: size(out{1,1}.max_risk_T,1)
        figure(f1)
        subplot(numel(idx2plot),1,i)
        post_val = out{te,1}.max_risk_T{idx2plot(i),path_idx ,1}.postNR_val;
        pre_val  = out{te,1}.max_risk_T{idx2plot(i),path_idx ,2}.preR_val;

        [h(te,i),p(te,i)] = kstest2(post_val,pre_val,'Tail','larger');

        cdfplot(post_val)
        hold
        cdfplot(pre_val)
        title(sprintf('%i p: %.4f  %s nCured: %i nImproved: %i',...
                        h(te,i),                                   ...
                        p(te,i),                                   ...
                        out{te,1}.risk_str_T{1}.Row{i},         ...
                       sum(~isnan(out{te,1}.max_risk_T{idx2plot(i),path_idx ,1}.postNR_val)),    ...
                       sum(~isnan(out{te,1}.max_risk_T{idx2plot(i),path_idx ,2}.preR_val))     ...
                      ),                                        ...
               'FontSize',14                                    ...
           );
       figure(f2)
       subplot(1,numel(idx2plot),i)
       g_idx = [zeros(size(pre_val)); ones(size(post_val)) ] ;
       %boxplot([pre_val;post_val],g_idx,'Labels',{'pre','post'})
       catNames = [repmat({'PRE resected'},size(pre_val),1); repmat({'POST'},size(post_val),1)] ;
       viol = violinplot([pre_val;post_val], catNames);
       viol(1).ViolinColor = [0 0 1];
       viol(2).ViolinColor = [1 0 0];
       
       c_bioName = (out{te,1}.risk_str_T{1}.Row{idx2plot(i)}); 
       title(c_bioName)
       MSize = 14;
       switch c_bioName
          case {'ARR'}
              ylabel(sprintf('max %s per subject',c_bioName))
          case {'PAC'}
              plot(1.5,0.37,'*','MarkerSize',MSize)
              ylabel(sprintf('max %s per subject',c_bioName))
          otherwise
              ylabel(sprintf('max strength per subject',c_bioName))
       end
    end   
    f2.PaperOrientation = 'landscape';
    set(f2, 'Position', get(0, 'Screensize'));
    print(sprintf('/home/matteo/Desktop/pics_4_result_section/actual_figures2/maxPreR_vs_Post_violinplot_%i',te),'-dpng')
    figure(f1)
    legend({'Post','Pre Resected'},'FontSize',18,'Location','northeast')
end



%% separability
cfg.nStep      = 100;
cfg.nStepQ     = 100;
cfg.normalize  = 1;
cfg.showpoints = 0;
cfg.spline     = 0;
%close all

nStep      = cfg.nStep;
nStepQ     = cfg.nStepQ;
normalize  = cfg.normalize;
showpoints = cfg.showpoints;



for te = 1 : size(out,1)
    figure
    for bio = 1 : size(out{te,1}.max_risk_T,1)

        post   = out{te,1}.max_risk_T{bio,path_idx ,1}.postNR_val;
        pre    = out{te,1}.max_risk_T{bio,path_idx ,2}.preR_val;
        maxBio = max([pre ; post]);
        minBio = min([pre ; post]);
        step   = (maxBio - minBio)/nStep; 

        binCenters = 0 : step : maxBio + 2*step;
      
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
            
            th_post          = max(post);
            [~,th_idx]       = min(abs(qpoints-th_post));
        
            subplot(4,2,bio)

            plot(qpoints,cs_post,'b')
            hold
            plot(qpoints,cs_pre,'r')
            %plot(qpoints,min(cs_post,cs_pre),'g--')

          

            if(showpoints)
                plot(binC,nC,'b*')
                plot(binI,nI,'ro')
            end
            
            th_y_val = 0;
            if(normalize)
                ylabel('subject density');
                th_y_val = 0.5;
            else
                ylabel('# subjects ');
                th_y_val = 10;
            end
            %line([th_post th_post],[0 th_y_val],'Color','green','LineStyle','--')
            line([qpoints(th_idx) qpoints(th_idx)],[0 th_y_val],'Color','green','LineStyle','--')
            
            overlap = cs_pre;
            overlap = overlap(th_idx:end);
            overlap = sum(overlap);
            
            
            
            
            title(sprintf('%s nCured: %i nImproved: %i \n >th: %i \n %i p: %.4f', ...
                          out{te,1}.risk_str_T{1}.Row{bio}, ...
                          sum(~isnan(post)),                ...
                          sum(~isnan(pre)),                 ...
                          sum(pre>th_post),                  ...
                          h(te,bio),                           ...
                          p(te,bio)                            ...
                         )                                  ...
                  )
    end
end
