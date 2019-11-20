function [idxArtefact ,idx_art_trial] = find_artefacts_epochs(s,labels,artefact_T)



% interval of the epoch
%s = outres.sampleinfo;


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
idxArtefact  = zeros(numel(labels),1);
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

       aux      = regexpi(labels,c_pattern);
       idxArtefact = ~cellfun(@isempty,aux);
       
        
end
