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
function subj_tbl = select_typeEPI(subj_tbl,typeEPI)


typeEPI_bin  = ~cellfun(@isempty,regexpi(subj_tbl.typeEPI,typeEPI));

idx_soi      = typeEPI_bin == 1;

subj_tbl = subj_tbl(idx_soi,:);
