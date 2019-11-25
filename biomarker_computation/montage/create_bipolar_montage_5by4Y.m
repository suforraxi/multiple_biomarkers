%create bipolar montage for grid 5x4 on the transversal direction 
% INPUT
% ch_label - channel names 
% data     - (ch X samples) data matrix 
% OUTPUT
% grid2use - grid labels available for which it was possible to compute the
%            bipolar montage
% outdata  - bipolar trasformation for the grid (i.e Gr1-Gr6) (channel X time sample)
function [grid2use,outdata]=create_bipolar_montage_5by4Y(ch_label,data)

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

  


    montage = [  1 6   ;...
                 6 11  ;...
                 11 16 ;...
                 2 7   ;...
                 7 12  ;...
                 12 17 ;...
                 3 8   ;...
                 8 13  ;...
                 13 18 ;...
                 4 9   ;...
                 9 14  ;...
                 14 19 ;...
                 5 10  ;...
                 10 15 ;...
                 15 20 ;...
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

