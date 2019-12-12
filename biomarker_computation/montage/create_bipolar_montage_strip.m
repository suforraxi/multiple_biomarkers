%create bipolar montage for strip 
% INPUT
% ch_label - channel names 
% data     - (ch X samples) data matrix 
% OUTPUT
% strip2use - strip labels available for which it was possible to compute
%             the bipolar montage
% outdata   - bipolar trasformation for the strip (i.e Str1-Str2) (channel X time sample)

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
function [strips2use,outdata]=create_bipolar_montage_strip(ch_label,data)


addhyphen = @(x) [x '-'];

strip2use = {};
outdata   = [];

% look for strip names, possible strip prefix are
strip_names  = {'Str','Rst','Riv','St','StrB','Sst','S1-'};


lookup4strip = regexpi(ch_label,'(Str|Rst|Riv|St|StrB|Sst|S1-)\d+');
strip_idx    = ~cell2mat(cellfun(@isempty,lookup4strip,'UniformOutput',false));
if(any(strip_idx))

    % extract only the strips
    strip_label  = ch_label(strip_idx);
    data         = data(strip_idx,:);

    % order the strips
    [~,reidx]    = sort(strip_label);
    strip_label  = strip_label(reidx);
    data         = data(reidx,:);

    % extract the prefix to construct the montage
    strip_prefix  = regexp(strip_label,'(?<stripname>\w+\D?[0]?)\d+','names');
    unique_prefix = {};

    for i = 1:numel(strip_prefix)
        unique_prefix{i} = strip_prefix{i}.stripname;
    end

    unique_prefix = unique(unique_prefix);
    
    %if(numel(unique_prefix)>1)
    %    error('more than one strip')
    %end
    
    for k = 1 : numel(unique_prefix)
    

        montage = { [1 2] ;...
                    [2 3] ;...
                    [3 4] ;...
                    [4 5] ;...
                    [5 6] ;...
                    [6 7] ;...
                    [7 8] ;...
                  };

        montage_M = zeros(size(montage,1),size(data,1));

            for i = 1:size(montage,1)

                c_mont      = montage{i};
                ch2check    = cell(size(c_mont));
                chAvailable = zeros(numel(strip_label),1);

                for j = 1:numel(c_mont)

                    ch2check{j} = sprintf( '%s%i' , unique_prefix{k}, c_mont(j) );

                    chAvailable = chAvailable | strcmp( strip_label ,ch2check{j});

                end

                curr_ch_idx = strcmp( strip_label , ch2check{1});

                if(sum(chAvailable) == numel(c_mont)) %all required channels are present

                   
                    chAvailable(curr_ch_idx) = 0;
                    montage_M(i,chAvailable) = -1;
                    montage_M(i,curr_ch_idx) = 1; 

                
                    strip2use{i,1} = strcat(ch2check{1},'-',ch2check{end});
                else
                    montage_M(i,:) = 0;
                    strip2use{i,1} = strcat(ch2check{1},'N-N',ch2check{end});
                end

            end


            montage2use = montage_M;
            
          
            outdata{k}    = montage2use * data;
            strips2use{k} = strip2use;
    end
else
    strips2use = {};
    outdata   = [];
end




