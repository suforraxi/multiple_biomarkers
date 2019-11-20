function [ res_channel, artefact_T] = get_metadata(dataDir,fileName)


[~,fileName,ext] = fileparts(fileName);
subFolder        = split(fileName,'_'); 

annotFile        = fullfile(dataDir,...
                     subFolder{1},...
                     subFolder{2},...
                     'ieeg',...
                      strcat(fileName,'_annotations.tsv'));


tsv_annots       = readtable(annotFile, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);
idx_resected     = strcmp(tsv_annots.type,'resected');
idx_edge         = strcmp(tsv_annots.type,'edge');

if(~isempty(idx_resected))
    res_channel      = tsv_annots.channel(idx_resected,1);
else
    res_channel      = [];
end
%edge_channel     = tsv_annots.channel(idx_edge,1);

idx_artefacts    = strcmp(tsv_annots.type,'artefact');
if(~isempty(idx_artefacts))
    artefact_T        = tsv_annots(idx_artefacts,:);
else
    artefact_T        = [];
end
