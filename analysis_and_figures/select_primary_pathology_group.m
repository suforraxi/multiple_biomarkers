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



function subj_tbl = select_primary_pathology_group(subj_tbl,group_id)

idx2keep = true(size(subj_tbl,1),1);
if (group_id ~= 0)
  
    idx2keep = subj_tbl.path_desc == group_id;
end

subj_tbl = subj_tbl(idx2keep,:);