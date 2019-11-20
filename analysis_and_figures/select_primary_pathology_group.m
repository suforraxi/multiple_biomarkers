function subj_tbl = select_primary_pathology_group(subj_tbl,group_id)

idx2keep = true(size(subj_tbl,1),1);
if (group_id ~= 0)
  
    idx2keep = subj_tbl.path_desc == group_id;
end

subj_tbl = subj_tbl(idx2keep,:);