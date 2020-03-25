function subject_above_th_per_biomarker()


load('/home/matteo/Desktop/mb_reviewers/avg_montage/figures/max_tables.mat')
subjAVG = aux_subject_above_th_per_biomarker(res_analysis)

load('/home/matteo/Desktop/mb_reviewers/figures/new_pool/max_tables.mat')
subjBip = aux_subject_above_th_per_biomarker(res_analysis)


comm_set = [];
diff_set = [];
for i = 1 : numel(subjAVG)

    comm_set{i} = intersect(subjAVG{i},subjBip{i});
    diff_set{i} = setdiff(subjAVG{i},subjBip{i});
end

 comm_set

function subj = aux_subject_above_th_per_biomarker(res_analysis)

idx_group    = 3;

idx_path     = 1;
idx_cured    = 1;
idx_improved = 2;

nBio = size(res_analysis.out{1}.max_T,1);
th   = zeros(nBio,1);


for i = 1 : nBio
    
    th(i)    = max(res_analysis.out{idx_group}.max_T{i,idx_path,idx_cured}.postNR_val);
    
    idx_subj = res_analysis.out{idx_group}.max_T{i,idx_path,idx_improved}.preR_val > th(i);
    
    subj{i} = res_analysis.out{idx_group}.max_T{i,idx_path,idx_improved}.subjName(idx_subj);
    
    
end