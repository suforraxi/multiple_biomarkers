% look if any channels of any trial contains artefacts
% INPUT
% s           - sample info information regarding the trials (what samples the trials was originated from the original recordings)            
% labels      - channel names considered 
% artefact_T  - table with the following variables
%               
%               type    name of the artefact (we have just generic artefact called 'artefact')        
%               start   sample time of the original data when the artefact starts        
%               stop    sample time of the original data when the artefact ends
%               channel channel name where the artefact was found
% OUTPUT
%
% idxArtefact   - array of booleans, the size of the array is equal to the
%                 number of channels. Every entry contains a boolean, 1
%                 artefact detected 0 no artefact
% idx_art_trial - array of booleans, the size of the array is equal to the
%                 number of trials. Every entry contains 1 if the
%                 corresponding trial has contains at least one channel
%                 with an artefact, 0 no artefact

%     Copyright (C) 2019 Matteo Demuru
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
function [idxArtefact ,idx_art_trial] = find_artefacts_epochs(s,labels,artefact_T)

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
