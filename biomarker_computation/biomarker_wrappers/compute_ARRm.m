%function to compute the ARRm values
%data - the recordings (channel X sample)
%windowL - number of sample in the windows 
% r1 - residual of the ARR model with order 1 (channel X #windows)
% r2 - residual of the ARR model with order 2 (channel X #windows)
% r3 - residual of the ARR model with order 3 (channel X #windows)
% r_clean - r3 but without possible artefacts (channel X # clean windows)
% h - DC offset (channel X #windows)
% an - 'slope of the DC' (channel X #windows)
% en - envelope (channel X #windows)
% avg_iPh - average instantaneous phase per windows (channel X #windows)
function [r1,r2,r3,r_clean,ARR,ARRm,h,an,en,avg_iPh]=compute_ARRm(data,windowL)
tic
h=0;an=0;en=0;avg_iPh=0;
%Add Stiliyan routines in order to use the ARR model
%addpath('/Users/matte/Desktop/git_rep/epi/matlab/stiliyan/matLabTools/')
% calculate ARR  
%windowL=40;
[~,r1,~] = get_residual(data,1,windowL);
[~,r2,~] = get_residual(data,2,windowL);
%it would be better to compute the DC (h) 'phase' (an) and envelope (en) on
% r_clean or clean the (h,an,en) later
%[~,r3,~,h,an,en,avg_iPh] = get_residual2(data,3,windowL);
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
toc
