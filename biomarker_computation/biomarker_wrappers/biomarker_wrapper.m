% wrapper for biomakers


function outdata = biomarker_wrapper(cfg,data)



if(cfg.cutTrials)
    
    cfgReTrials.trials  = cfg.trials_ct;
    cfgReTrials.length  = cfg.length_ct; %seconds of new trials
    cfgReTrials.overlap = cfg.overlap_ct;

    data = ft_redefinetrial(cfgReTrials,data);
end

% detrend and de mean
if(cfg.deT_deM)
    cfgPre.demean  = 'yes';
    cfgPre.detrend = 'yes';
    cfgPre.trial   = 'all';
    data = ft_preprocessing(cfgPre,data);
end

% notch filter
if(cfg.notch)

    cfgNotch.bsfilter = 'yes';
    cfgNotch.bsfreq   = cfg.notchBS;
    cfgNotch.trial    = 'all';
    
    data = ft_preprocessing(cfgNotch,data);
end

ntrial = numel(data.trial);




[~,outdata.hdr.datasetName,~] = fileparts(cfg.datasetName); 

outdata.label      = data.label;
outdata.time       = data.time;
outdata.fsample    = data.fsample;
outdata.sampleinfo = data.sampleinfo;


switch cfg.epiBio
    case 'ARR'
        [outdata.bio,outdata.extra] = wrapper_ARR(cfg,data);

    case 'PAC'
        [outdata.bio,outdata.extra] = wrapper_PAC(cfg,data);
    
    case 'PLI'    
        [outdata.bio,outdata.extra] = wrapper_PLI(cfg,data);
    
    case 'PLV'    
        [outdata.bio,outdata.extra] = wrapper_PLV(cfg,data);
    case 'H2'    
        [outdata.bio,outdata.extra] = wrapper_H2(cfg,data);
    case 'GC'    
        [outdata.bio,outdata.extra] = wrapper_GCtime(cfg,data);
    case 'sdDTF'
        [outdata.bio,outdata.extra] = wrapper_sdDTF(cfg,data);

end

outdata.type = cfg.epiBio; 
