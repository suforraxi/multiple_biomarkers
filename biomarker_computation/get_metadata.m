% function to obtain metadata regarding the recordings
% if a channel contains any artefact during the recordings
% 
% INPUT
% dataDir  - root folder for raw data (data should be organized in BIDS)
% fileName - filename for which to read the metadata (.vhdr) 
%
% OUTPUT
%
% artefact_T  - table with the following variables
%               
%               type    name of the artefact (we have just generic artefact called 'artefact')        
%               start   time of the original data when the artefact starts (seconds)        
%               stop    time of the original data when the artefact ends (seconds)
%               
%               channel channel name where the artefact was found

%     Copyright (C) 2020 Matteo Demuru
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

function [artefact_T] = get_metadata(dataDir,fileName)


[~,fileName,ext] = fileparts(fileName);
subFolder        = split(fileName,'_'); 

eventsFile        = fullfile(dataDir,...
                     subFolder{1},...
                     subFolder{2},...
                     'ieeg',...
                      strcat(fileName,'_events.tsv'));
channelsFile      = fullfile(dataDir,...
                     subFolder{1},...
                     subFolder{2},...
                     'ieeg',...
                      strcat(fileName,'_channels.tsv'));


tsv_channels     = readtable(channelsFile, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);              

sfreq            = tsv_channels.sampling_frequency(1);                  

tsv_annots       = readtable(eventsFile, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);

% convert from seconds to samples
tsv_annots.start = ceil(tsv_annots.start * sfreq);
tsv_annots.stop  = ceil(tsv_annots.stop * sfreq);

idx_artefacts    = strcmp(tsv_annots.type,'artefact');

if(~isempty(idx_artefacts))
    artefact_T        = tsv_annots(idx_artefacts,:);
else
    artefact_T        = [];
end
