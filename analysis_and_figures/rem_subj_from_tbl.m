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
function subj_T = rem_subj_from_tbl(subjName,subj_T)

for s = 1 : numel(subjName)
    
    idx = ~strcmp(subjName{s},subj_T.subjName);
    
    subj_T = subj_T(idx,:);
end