% select artefact free channels
%
% INPUT
% 
% subj_tbl   - table with the layout described in  create_summary_table_main
%              
% OUTPUT
%
% subj_tlb   - table with the layout described in
%              (create_summary_table_main) with only the entries with
%              channels without any artefacts
%            
function subj_tbl = select_artfree(subj_tbl)


idx_artfree = subj_tbl.artefact == 0;

subj_tbl = subj_tbl(idx_artfree,:);
