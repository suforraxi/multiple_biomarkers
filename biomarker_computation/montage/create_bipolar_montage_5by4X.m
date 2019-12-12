%create bipolar montage for grid 5x4 on the longitudinal direction 
% INPUT
% ch_label - channel names 
% data     - (ch X samples) data matrix 
% OUTPUT
% grid2use - grid labels available for which it was possible to compute the
%            bipolar montage
% outdata  - bipolar trasformation for the grid (i.e Gr1-Gr2) (channel X time sample)

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
function [grid2use,outdata]=create_bipolar_montage_5by4X(ch_label,data)

% extract grid channels
% look for grid different from 5x4

grid_idx = regexpi(ch_label,'Gr[0-9]?');
grid_idx = ~cellfun(@isempty,grid_idx);

if(sum(grid_idx)>20)
    error('grid bigger than 5x4')
end

if(any(grid_idx))
    grid_label = ch_label(grid_idx);
    data       = data(grid_idx,:)  ;

    % order the channels in ascending order and
    [~, reindex] = sort( str2double( regexp( grid_label, '\d+', 'match', 'once' )));
    grid_label   = grid_label(reindex);
    data         = data(reindex,:);

    

    montage = [  1 2   ;...
                 2 3   ;...
                 3 4   ;...
                 4 5   ;...
                 6 7   ;...
                 7 8   ;...
                 8 9   ;...
                 9 10  ;...
                 11 12 ;...
                 12 13 ;...
                 13 14 ;...
                 14 15 ;...
                 16 17 ;...
                 17 18 ;...
                 18 19 ;...
                 19 20 ;...
              ];

    montage2use = [];
    grid2use    = {};
    k           = 1; 
    for i = 1:size(montage,1)

        c_mont      = montage(i,:);
        ch2check    = cell(size(c_mont));
        chAvailable = zeros(numel(grid_label),1);

        for j = 1:numel(c_mont)

            if(any(cell2mat(regexpi(grid_label,'Gr0'))) && c_mont(j) < 10)

                prefix      = regexpi(grid_label,sprintf( 'Gr0%i' , c_mont(j) ),'match');
                idx_prefix  = ~cellfun(@isempty,prefix);

                if(all(~idx_prefix))
                    ch2check{j} = sprintf( 'Gr0%i' , c_mont(j));
                else

                    ch2check{j} = sprintf( '%s' , prefix{idx_prefix}{:});
                end
            else

                prefix      = regexpi(grid_label,'Gr','match');
                prefix      = prefix{1}{:};

                ch2check{j} = sprintf( '%s%i' , prefix,c_mont(j));
            end


            chAvailable = chAvailable | strcmp( grid_label ,ch2check{j});

        end

       if(sum(chAvailable) == numel(c_mont)) %all required channels are present

            curr_ch_idx1                        = strcmp( grid_label , ch2check{1}); 
            curr_ch_idxEnd                      = -1*strcmp( grid_label , ch2check{end}); 

            montage2use                         = [montage2use ; (curr_ch_idx1 + curr_ch_idxEnd )'] ;
            grid2use{i,1}                       = strcat(ch2check{1},'-',ch2check{end});
        else
            montage2use                         = [montage2use ; zeros(1,size(chAvailable,1))] ;
            grid2use{i,1}                       = strcat(ch2check{1},'N-N',ch2check{end});
       end


    end



    outdata    = montage2use * data;
else
    grid2use = []; 
    outdata  = [];
end

