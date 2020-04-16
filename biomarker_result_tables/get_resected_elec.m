% Extract the resected electrodes from *_eletrodes.tsv file
%
% 
% INPUT
% dataDir  - root folder for raw data (data should be organized in BIDS)
% fileName - filename for which to read the metadata (.vhdr) 
%
% OUTPUT
% resected_el - cell array with the name of the resected channels
%
%
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

function resected_el = get_resected_elec(dataDir,fileName)

[~,fileName,ext] = fileparts(fileName);
subFolder        = split(fileName,'_'); 

elecFile        = fullfile(dataDir,...
                     subFolder{1},...
                     subFolder{2},...
                     'ieeg',...
                      strcat(subFolder{1},'_',subFolder{2},'_ieeg_electrodes.tsv'));


tsv_t            = readtable(elecFile, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);
idx_resected     = strcmp(tsv_t.resected,'resected');
%idx_edge         = strcmp(tsv_annots.type,'edge');

if(~isempty(idx_resected))
    resected_el      = tsv_t.name(idx_resected,1);
else
    resected_el      = [];
end
