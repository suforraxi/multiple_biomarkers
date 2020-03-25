% Comparison of hippocampal electrodes and neocortical electrodes
%
% INPUT
% cfg   -  struct with the folloing fields
%             tbl2load               cell array with the filename of the table with
%                                    the results from the biomarker computation (see create_summary_table_main.m)
%             bioNames               cell array with the biomarker names  
%             alpha_level            double indicating the alpha level for the
%                                    Kolmogorov-Smirmov test used to compare hippocampal versus
%                                    neocortical biomakrker values
%             outFolderMaxComparison string with the folder name where tho save the result table
%
% OUTPUT 
% save a table with the biomarker name as row names the following variable
% 'pValue' which represents the p-value for the test 
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
function hippo_vs_neo_main(cfg)


tbl2load    = cfg.tbl2load;
bioNames    = cfg.bioNames; 
alpha_level = cfg.alpha_level;
outFolder   = cfg.outFolderMaxComparison;

h = zeros(numel(numel(tbl2load)),1);
p = zeros(numel(numel(tbl2load)),1);
for i = 1 : numel(tbl2load)

    % load table
    
    load(tbl2load{i})
    % use only temporal


    subj_tbl = select_typeEPI(subj_tbl,'T');

    % use only resected channels

    idx_resected = strcmp(subj_tbl.resected,'RES');

    subj_tbl     = subj_tbl(idx_resected,:);

    % look for neo-cortical and hippocampal channels
    idx_neo = [];
    idx_hip = [];
    switch cfg.montage
        case 'bipolar_two_directions'
            idx_neo = ~cellfun(@isempty,regexpi(subj_tbl.chName,'Gr.*-Gr.*'));

            idx_hip = ~cellfun(@isempty,regexpi(subj_tbl.chName,'(Str|Rst|Riv|St|StrB|Sst|S1-)[0]?(1|2|3)-(Str|Rst|Riv|St|StrB|Sst|S1-)[0]?(2|3|4)')) ;   
    
        case 'avg'
            idx_neo = ~cellfun(@isempty,regexpi(subj_tbl.chName,'Gr.*'));

            idx_hip = ~cellfun(@isempty,regexpi(subj_tbl.chName,'(Str|Rst|Riv|St|StrB|Sst|S1-)[0]?(1|2|3)')) ;   
    
            
    end
    
    neo_tbl = subj_tbl(idx_neo,:);

    hip_tbl = subj_tbl(idx_hip,:);

    neo_val = neo_tbl.biomarker;
    hip_val = hip_tbl.biomarker;


    [h(i),p(i,1)]   = kstest2(neo_val,hip_val,'Tail','larger');


    
end


p_tbl = table(num2str(p,'%1.5f'),'VariableNames',{'pValue'},'RowNames',bioNames');

writetable(p_tbl,fullfile(outFolder,'hippo_vs_neocort.tsv'),'FileType','text','Delimiter','tab','WriteRowNames',1);


