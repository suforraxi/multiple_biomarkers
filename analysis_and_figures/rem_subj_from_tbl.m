% remove subjects from table by name
function subj_T = rem_subj_from_tbl(subjName,subj_T)

for s = 1 : numel(subjName)
    
    idx = ~strcmp(subjName{s},subj_T.subjName);
    
    subj_T = subj_T(idx,:);
end