function pathstr_out = fullfile_cell(pathcell)
% function pathstr_out = fullfile_cell(pathcell)
%
% Does what MATLAB's fullfile() does, except with all path parts as items
% (chars) in a cell.
%
% Inputs:
% pathcell: cell containing all the parts of the path.
%
% Outputs:
% pathstr_out: path string.

for ii = 1:length(pathcell)
    if isempty(pathcell{ii}) % Convert empty elements to empty strings.
        pathcell{ii} = '';
    end
end

pathstr_out = pathcell{1};
if length(pathcell) == 1
    return
else
    for ii = 2:length(pathcell)
        pathstr_out = fullfile(pathstr_out, pathcell{ii});
    end
end
