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
