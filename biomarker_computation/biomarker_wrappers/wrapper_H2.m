% wrapper H2
function [ bio_vals, extra ] = wrapper_H2(cfg,data)

extra    = [];
bio_vals = [];

band_of_interest  = cfg.boi;

% zscore
for i = 1 : numel(data.trial)

    aux = data.trial{i};
    aux = zscore(aux,[],2);

    data.trial{i} = aux;

end
aux = [];

cfg_f            = [];
cfg_f.bpfilter   = 'yes';
cfg_f.bpfreq     = band_of_interest;
cfg_f.bpfilttype = 'fir';

f_data = ft_preprocessing(cfg_f, data);

[~,outdata.hdr.datasetName,~] = fileparts(cfg.datasetName); 

ntrial = size(f_data.trial,2);

for i = 1: ntrial
    
    m           = compute_h2(cfg,f_data.trial{i});
    aux.m       = m;
    aux.boi     = band_of_interest;
    extra{i}    = aux;
    
    bio_vals{i} = (sum(m,1) / (size(m,1)-1))'; % out strength
end



function m = compute_h2(cfg,a)
% a matrix (channels X time)
N     = size(a,1);
nch   = size(a,1);
m     = zeros(N); 


for i = 1 : nch
    for j = 1 : nch
       
       h2_aux = h2delay(a(i,:),a(j,:),cfg.T,cfg.n);
       m(i,j) =  max(h2_aux);
        
    end 
end


m(1:size(m,1)+1:end) = 0;
