% select type of epilepsy (temporal T / Extra Temporal / both)
function subj_tbl = select_typeEPI(subj_tbl,typeEPI)


typeEPI_bin  = ~cellfun(@isempty,regexpi(subj_tbl.typeEPI,typeEPI));

idx_soi      = typeEPI_bin == 1;

subj_tbl = subj_tbl(idx_soi,:);
