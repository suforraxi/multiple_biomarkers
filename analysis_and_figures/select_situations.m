function subj_tbl = select_situations(subj_tbl,reg_exp_str)


idx2keep  = ~cellfun(@isempty,regexp(subj_tbl.sitName,reg_exp_str));



subj_tbl = subj_tbl(idx2keep,:);