%% select artefact free channels
function subj_tbl = select_artfree(subj_tbl)


idx_artfree = subj_tbl.artefact == 0;

subj_tbl = subj_tbl(idx_artfree,:);
