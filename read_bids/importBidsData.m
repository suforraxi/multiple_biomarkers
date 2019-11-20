function [status,msg,cfg,data] = importBidsData(inDir)
try  
    status      = 0;
    msg         = [];
    cfg         = [];
    data        = [];
    
    hdrFile     = dir(fullfile(inDir,'*.vhdr'));
    hdrFile     = hdrFile(1).name;
    channelFile = dir(fullfile(inDir,'*channels.tsv'));
    channelFile = channelFile(1).name;
    annotFile   = dir(fullfile(inDir,'*annotations.tsv'));
    annotFile   = annotFile(1).name;

    
    cfg.datasetName = fullfile(inDir,hdrFile);
    cfg.channelFile = fullfile(inDir,channelFile);
    cfg.annotFile   = fullfile(inDir,annotFile);
    cfg.noArtefact  = 0;

    data = importSituation(cfg);
    
catch ME
   status = 1;
   msg    = sprintf('%s err:%s --func:%s',inDir,ME.message,ME.stack(1).name);
   data   = [];
   cfg    = [];
end