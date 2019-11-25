% compute sdDTF (see Mullen 2014 SIFT toolbox and  Zweiphenning 2019)
% 
% INPUT
% cfg - struct with the following field        
%           cfg.morder            :  model order to use for the estimation
%                                    of the multivariate model
%           cfg.freqBand          :  array specifying frequencies of interest
%           cfg.WindowLengthSec   :  window size for trial
%           cfg.WindowStepSizeSec :  step between windows
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
% bio_vals cell array with the out-strength computed on sdDTF matrix 
%          (average over functional connectivity rows extra.m) computed for every channel and every trial
%
% extra    struct with the following fields
%               extra.m      directed connectivity matrix 
%               extra.morder model order used to estimate the model
%               extra.freq   frequency where the model was estimated
%               extra.IC     model fitting information (see SIFT toolbox) 
%               extra.CONN   connectivity information  (see SIFT toolbox)
% 

function [ bio_vals, extra ] = wrapper_sdDTF(cfg,data)

extra    = [];
bio_vals = [];


ntrial = size(data.trial,2);

for i = 1: ntrial
    try
        hdr         = data.hdr;
        
        hdr.nSamples = size(data.trial{i},2);
        hdr.nTrials  = 1;
        % zscore
        zData       = zscore(data.trial{i},[],2);
        res         = compute_sdDTF(cfg,hdr,zData);

        aux.m      = res.m;
        aux.morder = res.morder;
        aux.freq   = res.freq;
        aux.IC     = res.IC;
        aux.CONN   = res.CONN;
        
        extra{i}    = res;

        bio_vals{i} = (nansum(aux.m,1) / (size(aux.m,1)-1))'; % out strength
    catch ME
        aux.m      = [];
        aux.morder = [];
        aux.freq   = [];
        aux.IC     = [];
        aux.Conn   = [];
        extra{i}   = [];

        bio_vals{i} = []; 
        
    end
end

% code adapted from code used by in the manuscript Zweiphenning 2019 for
% the sdDTF
%
function res = compute_sdDTF(cfg,hdr,data)

res = [];
% transform the data into eeglab data structure
EEG = fieldtrip2eeglab(hdr,data);
                
VERBOSITY_LEVEL = 2; 


%prep_eeg =pre_prepData(EEG,'VerbosityLevel',VERBOSITY_LEVEL,'SignalType',{'Channels'},'Detrend',[],  ...
%                        'NormalizeData',[]);

EEG = pop_pre_prepData(EEG,'nogui',...
                            'VerbosityLevel',VERBOSITY_LEVEL,...
                            'SignalType',{'Channels'},...
                            'Detrend',1,...
                            'resetConfigs',true,...
                            'badsegments',[],...
                            'newtrials',[],...
                            'NormalizeData',{'Method','time'} ...
                            );

% model estimation

%EpochTimeRange      = [0 5];                 % this is the time range (in seconds) to analyze (relative to event at t=0)
WindowLengthSec     = cfg.WindowLengthSec;%[5];        % sliding window length in seconds
WindowStepSizeSec   = cfg.WindowStepSizeSec;%[5];      % sliding window step size in seconds
SelAlg              = {'ARfit'};      %,'Group Lasso (DAL/SCSA)','Group Lasso (ADMM)'};        % selection algorithms to determine model order
%NewSamplingRate     = [];                           % new sampling rate (if downsampling)
GUI_MODE            = 'nogui';                      % whether or not to show the Graphical User Interfaces. Can be 'nogui' or anything else (to show the gui)
VERBOSITY_LEVEL     = 2;                            % Verbosity Level (0=no/minimal output, 2=graphical output)

morder = cfg.morder;
freq   = cfg.freqBand;

[EEG.CAT.IC, EEG.CAT.configs.est_selModelOrder]  = est_selModelOrder(EEG,...
    'modelingApproach',{'Segmentation VAR','algorithm',SelAlg(1),'winlen',WindowLengthSec,'winstep',WindowStepSizeSec,...
                    'prctWinToSample',100,'normalize',[],'detrend',[],'setArgDirectMode',true,'verb',VERBOSITY_LEVEL,...
                    'ModelOrder',morder},...
                    ...
                    'morderRange',[morder morder],...
                    'icselector',{'sbc','aic', 'hq','fpe'},...
                    'downdate',1,...
                    'runPll',{'ProfileName','local','NumWorkers',6},...
                    'plot',[],...
                    'verb',VERBOSITY_LEVEL);


EEG = pop_est_fitMVAR(EEG,GUI_MODE,EEG.CAT.IC.modelFitting.modelingArguments,'ModelOrder',morder);

EEG = pop_est_mvarConnectivity(EEG,GUI_MODE, ...
                            'connmethods',{'DTF' 'nDTF' 'dDTF' 'dDTF08' 'ffDTF'}, 'absvalsq',true,'spectraldecibels',false,   ...
                            'freqs',freq ,'verb',VERBOSITY_LEVEL);
                        
                        
m = EEG.CAT.Conn.DTF;
m = m./(sum(sum(sum(m)))); % divide by total connectivity
% mean over the frequency band
m = sum(m,3);               
m(1:(size(m,1)+1):end) = 0; % diagonal to zero
m = m./(size(m,1)-1);

res.m      = m;
res.morder = morder;
res.freq   = freq;
res.IC     = EEG.CAT.IC.modelFitting.modelingArguments;
res.CONN   = EEG.CAT.Conn;
                        