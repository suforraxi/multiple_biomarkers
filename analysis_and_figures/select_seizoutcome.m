% select the subject according to a regular expression 
% on the type of seizure outcome i.e. '1a\w*'

% INPUT
% 
% subj_tbl   - table with the layout described in  (create_summary_table_main)
% sf_var     - variable in the subj_tbl to consider for the selection
% regexp_STR - regular expression to describe which patients will be kept
%              in the table 
% OUTPUT
%
% subj_tlb   - table with the layout described in
%              (create_summary_table_main) with only the entries which satisfy the regular expression 

function subj_tbl = select_seizoutcome(subj_tbl,sf_var,regexp_STR)

seizout_bin = ~cellfun(@isempty,regexpi(eval(sprintf('subj_tbl.%s',sf_var)),regexp_STR));

idx_so      = seizout_bin == 1;

subj_tbl = subj_tbl(idx_so,:);