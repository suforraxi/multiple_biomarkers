% wrapper for ARR
function [ bio_vals, extra ] = wrapper_ARR(cfg,data)

bio_vals = [];
extra    = [];
ntrial = size(data.trial,2);

% zscore
for i = 1 : numel(data.trial)

    aux = data.trial{i};
    aux = zscore(aux,[],2);

    data.trial{i} = aux;

end
aux = [];


for i = 1: ntrial
    
    
    
    [~,...
     ~,...
     aux.r3,...
     aux.rclean,...
     aux.ARR,...
     aux.ARRm,...
     ~,~,~,~ ...
     ]...
     = compute_ARRm(data.trial{i},cfg.windowL);

    bio_vals{i} = aux.ARRm';
    
    aux.r3      = aux.r3';
    aux.rclean  = aux.rclean';
    aux.ARR     = aux.ARR';
    aux.ARRm    = aux.ARRm';
    
    extra{i}    = aux;
end

