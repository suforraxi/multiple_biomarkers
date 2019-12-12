% wrapper H2 (Kalitzin 2006)
%
% INPUT
% cfg - struct with the following field        
%        cfg.boi : [x y] frequency band boundaries to filter the signals before to compute H2
%        cfg.T   : specify an array of delays to compute H2 (see external/h2delay.m)
%        cfg.n   : number of bins see (external/h2delay.m)
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
% bio_vals cell array with the H2 out-strength (average over functional connectivity rows extra.m) 
%          computed for every channel and every trial
%
% extra    struct with the following fields
%           extra.m    functional connectivity matrix of H2 values (maximum value across delays)   
%           extra.boi  frequency band boundaries used to filter the signals before to compute H2 
% 

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
       m(i,j) =  max(h2_aux); %save maximum value across delays
        
    end 
end


m(1:size(m,1)+1:end) = 0;
