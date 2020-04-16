% generic wrapper for biomakers 
% all the common pre-processing steps are done in this function like
% 1) cut the last minute
% 2) demean and detrend
% 3) apply notch filter

% INPUT
% cfg struct with the following fields 
%         datasetName  : file name (.vhdr) that will be analysed (see BIDS format)
%         channelFile  : channel file name relative to the channel file in BIDS (see BIDS format) 
%         annotFile    : annotation file we used a custom annotation file for ioECoG (see BIDS + https://github.com/suforraxi/ieeg_respect_bids )
%         noArtefact   : field to remove it is not used anymore
%         inDir_data   : input root directory for the data in BIDS
%         subj_info_F  : file name containing a subject information table with variables specified in (batch_compute_different_biomarkers.m)
%         cutLast      : 1 in order to cut the length amount of data (last minute in our case reduce propofol effect) of the recordings 
%                        0 otherwise    
%         trials       : fieldtrip field used to implement the cut (see ft_preprocessing) 
%         length       : fieldtrip field used to implement the cut (see ft_preprocessing) 
%         overlap      : fieldtrip field used to implement the cut (see ft_preprocessing) 
%         cutTrials    : 1 in order to redefine the trials length 
%                        0 otherwise
%         trials_ct    : fieldtrip field used to implement the cut (see ft_redefinetrials) 
%         length_ct    : fieldtrip field used to implement the cut (see ft_redefinetrials) 
%         overlap_ct   : fieldtrip field used to implement the cut (see ft_redefinetrials) 
%         deT_deM      : 1 to apply detred and demean (see ft_preprocesing)
%                        0 otherwise               
%         notch        : 1 to apply notch filter
%                        0 otherwise   
%         notchBS      : frequency interval where to apply the notch filter 
%                        [low high] (see ft_preprocessing)
%         montage      : 'bipolar_two_directions' apply the bipolar montage for
%                         grid 5x4 and strip 1x6 or 1x8 (custom function to compute montage see /montage/ folder)
%         outdir_combi : folder name where to save the results (outres see below)
%         errorFile    : file name where to save the failures
%         epiBio       : biomaker name to compute it could be (ARR / PAC / PLV / PLI / H2 / GC / sdDTF )
% 
%         extra field required, depending on the biomarker
%          ARR  (see wrapper ARR.m)
%                   windowL : length in samples of the slinding window used to compute the ARR   
%          PAC  (see wrapper_PAC.m)    
%                   lb: low frequency band boundaries [x y] used to estimate the phase   
%                   hb: high frequency band boundaries [x y] used to estimate the
%                       amplitude envelope
%          PLI  (see wrapper_PLI.m)
%                   boi : [x y] frequency band boundaries to filter the signals before to compute PLI
%          PLV  (see wrapper_PLV.m)
%                   boi : [x y] frequency band boundaries to filter the signals before to compute PLV
%          H2   (see wrapper_H2.m) 
%                   boi : [x y] frequency band boundaries to filter the signals before to compute H2
%                   T   : specify an array of delays to compute H2 (see external/h2delay.m)
%                   n   : number of bins see (external/h2delay.m)
%          GC   (see wrapper_GCtime.m)
%                   momax : maximal model order to try
%          sdDTF (see wrapper_sdDTF.m)
%                   morder            :  model order 
%                   freqBand          :  array specifying frequencies of interest
%                   WindowLengthSec   :  window size for trial
%                   WindowStepSizeSec :  step between windows
% 
%
% data  - fieldtrip data structure (see http://www.fieldtriptoolbox.org/)
%          
%             data.trial
%             data.time
%             data.fsample
%             data.label
%             data.sampleinfo
%               
%
% OUTPUT
%
% outdata - The generic result is saved in a matlab structure with the following fields
% 
%                           outdata.hdr          - containing the datasetName from which the
%                                                 biomarker is computed
%                           outdata.label        - cell array containing the channel names for which
%                                                 the biomarker is computed 
%                           outdata.time         - cell array containing the time for each trial in
%                                                 which the original data is divided
%                           outdata.fsample      - sample frequency of the data
%                           outdata.sampleinfo   - sample information corresponding to the samples
%                                                 in the original recordings from which each trial
%                                                 computed
%                           outdata.bio          - cell array with an entry per trial of the input
%                                                 data. Each entry corresponds another cell array
%                                                 with the computation of the biomarker per each
%                                                 channel specified by label
%                                                 (or biomaker statistic, like the strenght for bivariate/ multivariate methods)
%                           outdata.extra        - struct with extra results depending on the
%                                                 biomaker (see biomarker_wrapper specific function)
%                           outdata.type         - name of the biomaker computed
           

%     Copyright (C) 2019 Matteo Demuru
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
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




[~,fName,~]                    = fileparts(cfg.datasetName);
artefact_T                     = get_metadata(cfg.inDir_data,fName);
[idxChArtefact ,idx_art_trial] = find_artefacts_epochs(data.sampleinfo,data.label,artefact_T);

cfgCH.channel = data.label(~idxChArtefact);
data          = ft_preprocessing(cfgCH,data);  

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
