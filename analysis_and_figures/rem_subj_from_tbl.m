% remove subjects from table by name

% INPUT
% 
% subjName   - cell array with the names ('RESPXXXX') of subjects to remove
%              from the subj_tbl
% subj_tbl   - table with the layout described in  create_summary_table_main
% 
%               
% OUTPUT
%
% subj_tlb   - table with the layout described in
%              (create_summary_table_main) without the subjects specified
%              in subjName

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
function subj_T = rem_subj_from_tbl(subjName,subj_T)

for s = 1 : numel(subjName)
    
    idx = ~strcmp(subjName{s},subj_T.subjName);
    
    subj_T = subj_T(idx,:);
end