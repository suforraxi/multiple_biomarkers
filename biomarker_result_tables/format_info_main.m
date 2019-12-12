% count how many grid 5 by 4 +  zero strip / one strip / two strips or more
% 
% INPUT
% cfg struct with the following fields
%  bidsDir  - folder name of the raw data organized in BIDS format (see https://github.com/suforraxi/ieeg_respect_bids)
%  InfoF    - file name of the information table (see batch_compute_different_biomarkers)
% 
% OUTPUT
% sit_info_T - table with the following variables 
%                          nameXsit      - name of the situation file considered                     
%                          subjNameXsit  - coded name of the subject (RESPXXXX)              
%                          formatXsit    - format of the grid/strip used for the recording (see create_summary_table_main for details)            
%                          nGridXsit     - number of grid used in the situation 
%                          nStripXsit    - number of strip used in the situation
%                          nExcXsit      - number of format exceptions

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
function sit_info_T = format_info_main(cfg)

bidsDir  = cfg.bidsDir; 
% load information file with subjects
InfoF    = cfg.InfoF; 

info_T   = readtable(InfoF,'Delimiter','\t','FileType','text','ReadRowNames',1,'ReadVariableNames',1);

%load(tleInfoF);

subjName     = info_T.Properties.RowNames;
nameXsit     = [];
subjNameXsit = [];
nGridXsit    = [];
nStripXsit   = [];
nExcXsit     = [];
formatXsit   = [];

sit_c        = 0;
for i = 1 : numel(subjName)
    % find situations for subject
    subjDir  = fullfile(bidsDir,strcat('sub-',subjName{i}));
    sitList = dir(fullfile(subjDir,'*SITUATION*'));
    
    for j = 1 : numel(sitList)
        
        sit_c   = sit_c + 1;
        
        sitName = sitList(j).name;
        
        
        subjNameXsit{sit_c,1} = subjName{i};
        
        % find the json with the Format information
        vhdrF = dir(fullfile(subjDir,sitName,'ieeg','*vhdr'));
        if(~isempty(vhdrF))
            
            jsonF = replace(vhdrF(1).name,'.vhdr','.json');
            jsonF = fullfile(vhdrF(1).folder,jsonF);
            
            if(isfile(jsonF))
                
                [~,nameXsit{sit_c,1},~] = fileparts(jsonF);
                
                jsonSideCar = loadjson(jsonF);
                format_str  = jsonSideCar.iEEGElectrodeGroups;
                % Format;Gr[4x5];(Str|Rst|Riv|St|StrB|Sst|S1-)[];(Str|Rst|Riv|St|StrB|Sst|S1-)[];*
                % Format;Exception
                
                format_split = regexpi(format_str,';','split');
                format_split = format_split(2:end);
                      
                % check if contains grid 4 by 5
                look_gr    = regexpi(format_split,'Gr\[4x5\]','match');
                hasGrid    = ~cellfun(@isempty,look_gr);
                nGrids     = sum(hasGrid); 
                % check if contains strips
                look_strip = regexpi(format_split,'(Str|Rst|Riv|St|StrB|Sst|S1-)\[.*\]','match');
                hasStrips  = ~cellfun(@isempty,look_strip);
                nStrips    = sum(hasStrips);
                % check if contains exception
                look_exc   = regexpi(format_split,'Exception','match');
                hasExc     = ~cellfun(@isempty,look_exc);
                nExc       = sum(hasExc);

                nGridXsit(sit_c,1)  = nGrids;
                nStripXsit(sit_c,1) = nStrips;
                nExcXsit(sit_c,1)   = nExc; 
                formatXsit{sit_c,1} = format_str;  
            end
        end
    end
end

sit_info_T = table(nameXsit,subjNameXsit,formatXsit,nGridXsit,nStripXsit,nExcXsit);

