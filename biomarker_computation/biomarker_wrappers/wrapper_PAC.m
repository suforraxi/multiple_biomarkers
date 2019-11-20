% wrapper PAC 
function [ bio_vals, extra ] = wrapper_PAC(cfg,data)

extra    = [];
bio_vals = [];


% zscore
for i = 1 : numel(data.trial)

    aux = data.trial{i};
    aux = zscore(aux,[],2);

    data.trial{i} = aux;

end
aux = [];



low_band  = cfg.lb;
high_band = cfg.hb;

cfg_LB            = [];
cfg_LB.bpfilter   = 'yes';
cfg_LB.bpfreq     = low_band;
cfg_LB.bpfilttype = 'fir';

low_data = ft_preprocessing(cfg_LB, data);

cfg_HB            = [];
cfg_HB.bpfilter   = 'yes';
cfg_HB.bpfreq     = high_band;
cfg_HB.bpfilttype = 'fir';

high_data = ft_preprocessing(cfg_HB, data);

[~,outdata.hdr.datasetName,~] = fileparts(cfg.datasetName); 

ntrial = size(data.trial,2);

for i = 1: ntrial

    hA = hilbert(high_data.trial{i}');
    Lp = hilbert(low_data.trial{i}');

    Lp = Lp./abs(Lp);
    
    a = diag(Lp'*abs(hA));
    b = diag(abs(Lp')*abs(hA));

    bio_vals{i}  = abs(a./b); 
    pac          = a./b;
    
    extra_s.pac = pac;
    extra_s.lf  = low_band;
    extra_s.hf  = high_band;
    
    extra{i}    = extra_s;
     
end



