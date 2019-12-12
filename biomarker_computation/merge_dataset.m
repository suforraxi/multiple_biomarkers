% merge the data fields of data1 and data2 in merge_data
% used to merge grid data in different directions and strip data
% INPUT
% data1 and data2 - fieldtrip data struct to merge
% OUTPUT
% merged_data     - fieldtrip struct containing the merged fields
%
% Note it removes channels not present indicated by N-N in the name for
% bipolar montage

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
function merged_data = merge_dataset(data1,data2)

    merged_data = [];

    if(isempty(data1) && ~isempty(data2))
        merged_data = data2;
    end
    
    if(isempty(data2) && ~isempty(data1))
        merged_data = data1;
    end
    
    if(~isempty(data1) && ~isempty(data2))
        
        merged_data.sampleinfo = data1.sampleinfo;          
        
        merged_data.label   = [data1.label ; data2.label];
        merged_data.time    =  data1.time;

        for i = 1 : numel(data1.trial)
            merged_data.trial{i}   = [data1.trial{i} ; data2.trial{i} ];
        end

        merged_data.fsample = data1.fsample;
    end
    if(~isempty(merged_data))
        idx2rem = [];
        for i = 1 : numel(merged_data.label)

            if(contains(merged_data.label{i},'N-N'))

                idx2rem = [idx2rem i];
            end
        end

        merged_data.label(idx2rem) = [];
        for i = 1 : numel(merged_data.trial)

            merged_data.trial{i}(idx2rem,:)= []; 

        end
    
    end