% select the subject according to a regular expression 
% on the type of seizure outcome i.e. '1a\w*'
function subj_tbl = select_seizoutcome(subj_tbl,sf_var,regexp_STR)

seizout_bin = ~cellfun(@isempty,regexpi(eval(sprintf('subj_tbl.%s',sf_var)),regexp_STR));

idx_so      = seizout_bin == 1;

subj_tbl = subj_tbl(idx_so,:);