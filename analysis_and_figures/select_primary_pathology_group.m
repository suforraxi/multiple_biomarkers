% select primary pathology subjects (see create_summay_table_main)
%
% INPUT
% 
% subj_tbl   - table with the layout described in  create_summary_table_main
% group_id   - integer corresponding to the primary pathology index to keep (see create_summary_table_main) 
%               
% OUTPUT
%
% subj_tlb   - table with the layout described in
%              (create_summary_table_main) with only the entries with the
%              primary pathology specified by group_id

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

function subj_tbl = select_primary_pathology_group(subj_tbl,group_id)

idx2keep = true(size(subj_tbl,1),1);
if (group_id ~= 0)
  
    idx2keep = subj_tbl.path_desc == group_id;
end

subj_tbl = subj_tbl(idx2keep,:);