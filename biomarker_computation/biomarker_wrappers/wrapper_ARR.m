% wrapper for ARR (Geertsema 2017)
%
% INPUT
% cfg - struct with the following field
%        cfg.windowL : window length to compute compute the residual expressed in number of sample 
% 
% data - fieldtrip data structure
%         data.trial
%         data.time
%         data.fsample
%         data.label
%         data.sampleinfo
%
%
%
% OUTPUT
% bio_vals cell array with the ARRm computed for every channel and every
%          trial
% extra    struct with the following fields
%           r1      - residual of the ARR model with order 1 (channel X #windows)
%           r2      - residual of the ARR model with order 2 (channel X #windows)
%           r3      - residual of the ARR model with order 3 (channel X #windows)
%           r_clean - r3 but without possible artefacts (channel X # clean windows)
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
     aux.ARRm...
     ]...
     = compute_ARRm(data.trial{i},cfg.windowL);

    bio_vals{i} = aux.ARRm';
    
    aux.r3      = aux.r3';
    aux.rclean  = aux.rclean';
    aux.ARR     = aux.ARR';
    aux.ARRm    = aux.ARRm';
    
    extra{i}    = aux;
end

%function to compute the ARRm values
% data    - the recordings (channel X sample)
% windowL - window length to compute compute the residual expressed in number of sample 
%
% r1      - residual of the ARR model with order 1 (channel X #windows)
% r2      - residual of the ARR model with order 2 (channel X #windows)
% r3      - residual of the ARR model with order 3 (channel X #windows)
% r_clean - r3 but without possible artefacts (channel X # clean windows)

function [r1,r2,r3,r_clean,ARR,ARRm]=compute_ARRm(data,windowL)


% see /external/get_residual.m
[~,r1,~] = get_residual(data,1,windowL);
[~,r2,~] = get_residual(data,2,windowL);

[~,r3,~] = get_residual(data,3,windowL);
r_clean=zeros(size(r3));

%clean the ARR see Geertsema 2017 Clin. Neurophysiology
nchs=size(data,1);
for c=1:nchs
    % residual decline over residual from order 1 and 2 models
    rn = [r1(c,:);r2(c,:)];
    D = -2*diff(rn)./(rn(1,:)+rn(2,:));

    % remove high r3 values probably resulting from artefacts
    r_clean(c,:) = r3(c,:);
    p95 = prctile(r_clean(c,:),95);
    w=length(D);
    while w>0
        
        if r_clean(c,w)>p95 && D(w)<0.9
            %remove possible artefacts at w
            %and two neighbors windows
            %checking for border conditions
            
            r_clean(c,w) =NaN;%NaN artecfact
            
            if(w~=length(D))
                r_clean(c,w+1) =NaN;%NaN artecfact
            end
            
            if(w~=1)
                r_clean(c,w-1) =NaN;%NaN artecfact
            end
            
            w=w-1;
        
        end
        w=w-1;
    end
    
    if(any(data(c,:),2))
        ARR(c) = std(r3(c,:))/mean(r3(c,:));
        ARRm(c) = nanstd(r_clean(c,:))/nanmean(r_clean(c,:));
    else
        ARR(c) = NaN;
        ARRm(c) = NaN;
    end
    
    
end

