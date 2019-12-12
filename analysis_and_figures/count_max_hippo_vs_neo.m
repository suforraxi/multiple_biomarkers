% Count how many times the maximum is on the grid (neocortical channel) or on the
% first three channels of the strip (hippocampal channel)
%
%
% INPUT 
% res_analysis   struct with the following fields
%                
%                out:               cell array containing struct for per type of epilepsy / seizure outcome group 
%                                   Each struct as a field max_T storing a 3-dimension cell array  
%                                   organized as biomaker X pathology group X seizure outcome group
%                                   each cell contains a table with the following variables
% 
%                                       subjName     - coded subject name
%                                       postNR_sit   - situation name from post-resection where
%                                                      the maximum value of the biomarker is found 
%                                       post_chName  - channel name where the the maximum value of the biomarker is found
%                                       postNR_val   - biomarker value corresponding to the maximum post-resection 
%                                       preR_sit     - situation name from pre-resection recordings where
%                                                      the maximum value of the biomarker is found across resected channels 
%                                       preR_chName  - channel name where the the maximum value
%                                                      of the biomarker is found across resected channels 
%                                       preR_val     - biomarker value corresponding to the
%                                                      maximum across pre-resection resected channels
%                                       preNR_sit    - situation name from pre-resection recordings where
%                                                      the maximum value of the biomarker is found across not resected channels
%                                       preNR_chName - channel name where the the maximum value
%                                                      of the biomarker is found across not resected channels
%                                       preNR_val    - biomarker value corresponding to the
%                                                      maximum across pre-resection not resected channels
% 
% 
%                seizOUT:          cell array with strings describing the seizure
%                                  outcome group contained in out i.e.  {'1a_AED_stop\w*', '1(a|b)\w*'} for cured patients and improved patients
% 
%                typeEPI:          cell array with strings describing the type of epilepsy group i.e. {'Joint_T_E'  'Temporal'  'Extra_Temporal'}  
%                path_group_label: cell array with a string describing the pathology group  i.e. {'all_primary_pathologies'}
%             
%                seizOutVariable:  cell array with a string describing the seizure outcome variable consideres i.e. {'description_sf_1y'}
%
%
% outFolder      string indicating the folder where the output will be saved
% bioNames       cell array with the name of the biomarkers considered i.e.
%                {'ARR','PAC','PLV','PLI','H2','GC','sdDTF'}
%                NB biomarkers should be ordered as they are store in out.max_T 
%
%
% OUTPUT
%
% save table (r_tbl) with the biomarker name as row names and two variables
%                     
%                      MaxOnHippocampus     - number of subjects for whom the maximum across all channels/situations during the pre-resection recording was found in the hippocampal channels
%                      subj_greater_than_th - total number of subjects for whom the maximum was higher than the threshold computed as biomarker reference using post-resection recordings
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
function count_max_hippo_vs_neo(res_analysis,outFolder,bioNames)



% temporal subjects

temporal_idx = 2;
nBiomarkers  = size(res_analysis.out{temporal_idx}.max_T,1);
path_idx     = 1;
improved_idx = 2;
cured_idx    = 1; 

totNhit  = zeros(nBiomarkers,1);
counts   = zeros(nBiomarkers,1);

for b = 1 : nBiomarkers  
    
   maxCured   =  max(res_analysis.out{temporal_idx}.max_T{b,path_idx,cured_idx}.postNR_val);
   
   preR_val   = res_analysis.out{temporal_idx}.max_T{b,path_idx,improved_idx}.preR_val;
   
   idx_hit    = preR_val > maxCured;
   
   totNhit(b) = sum(idx_hit);
   
   chName     = res_analysis.out{temporal_idx}.max_T{b,path_idx,improved_idx}.preR_chName(idx_hit);

   idx_hip    = ~cellfun(@isempty,regexpi(chName,'(Str|Rst|Riv|St|StrB|Sst|S1-)[0]?(1|2|3)-(Str|Rst|Riv|St|StrB|Sst|S1-)[0]?(2|3|4)')) ; 
   
   counts(b)  = sum(idx_hip);
end

r_tbl = table(counts,totNhit,'VariableNames',{'MaxOnHippocampus','subj_greater_than_th'},'RowNames',bioNames');

writetable(r_tbl,fullfile(outFolder,'max_hippo_vs_neocort.tsv'),'FileType','text','Delimiter','tab','WriteRowNames',1);