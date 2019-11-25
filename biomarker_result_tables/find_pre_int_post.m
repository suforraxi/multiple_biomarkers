% find pre / intermediate / post situations  
%
% INPUT
% situations - cell array containing situations names
% 
% 
% OUTPUT
%
% pre    - cell array with the name of the pre-resection situations
% inter  - cell array with the name of the intermidiate situations
% post   - cell array with the name of the post-resection situations
function [pre,inter,post] = find_pre_int_post(situations)

pre   = [];
inter = [];
post  = [];

if (~isempty(situations))
    sitName   = regexpi(situations,'\D*(?<num>\d)(?<letter>\D)(?<rest>\w*)','names');
    sitPhases = zeros(1,numel(sitName));

    for i = 1 : numel(sitName)

        sitPhases(i) = str2num(sitName{i}.num);
    end

    lastSit = max(sitPhases);
    if(lastSit==1)
        lastSit = 500; % to fix in a better way
    end
    
    pre   = contains(situations,'SITUATION1');
    
    post  = contains(situations,sprintf('SITUATION%i',lastSit));
    
    inter = ~(pre | post);

end
