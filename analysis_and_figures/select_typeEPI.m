% select type of epilepsy (temporal T / Extra Temporal / both)
%
% INPUT
% 
% subj_tbl   - table with the layout described in  (create_summary_table_main)
% typeEPI    - regular expression to describe which patients will be kept
%              in the table 
% OUTPUT
%
% subj_tlb   - table with the layout described in
%              (create_summary_table_main) with only the entries which satisfy the regular expression

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
function subj_tbl = select_typeEPI(subj_tbl,typeEPI)


typeEPI_bin  = ~cellfun(@isempty,regexpi(subj_tbl.typeEPI,typeEPI));

idx_soi      = typeEPI_bin == 1;

subj_tbl = subj_tbl(idx_soi,:);
