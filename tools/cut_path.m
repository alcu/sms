function path_out = cut_path(path_in)
% function path_out = cut_path(path_in)
%
% For a given path_in string, outputs just the directories and files in the
% path.
%
% Inputs:
% path_in: string portraying some path. Make sure you don't have any carriage
% returns at the end of path_in.
%
% Outputs:
% path_out: cell, where each element of the cell is a directory or file, in
% order of path_in.

idx = find(path_in == filesep);

nfseps = length(idx);

if nfseps == 0
    path_out{1} = path_in;
    return
end

if idx(1) == 1 & idx(end) == length(path_in) % Beg and last chars of path_in are fileseps.
    path_out = cell(1, nfseps - 1);
    for ii = 1:length(path_out)
        path_out{ii} = path_in(idx(ii) + 1 : idx(ii + 1) - 1);
    end
elseif idx(1) == 1 % Only beginning is a filesep.
    path_out = cell(1, nfseps);
    for ii = 1:length(path_out) - 1
        path_out{ii} = path_in(idx(ii) + 1 : idx(ii + 1) - 1);
    end
    path_out{end} = path_in(idx(end) + 1:end);
elseif idx(end) == length(path_in) % Only end is a filesep.
    path_out = cell(1, nfseps);
    path_out{1} = path_in(1:idx(1) - 1);
    for ii = 2:length(path_out)
        path_out{ii} = path_in(idx(ii - 1) + 1 : idx(ii) - 1);
    end
else % Beg and end are not fileseps.
    path_out = cell(1, nfseps + 1);
    path_out{1} = path_in(1:idx(1) - 1);
    for ii = 2:length(path_out) - 1
        path_out{ii} = path_in(idx(ii - 1) + 1 : idx(ii) - 1);
    end
    path_out{end} = path_in(idx(end) + 1:end);
end

empty_items = [];
for ii = 1:length(path_out)
    if isempty(path_out{ii})
        empty_items = [empty_items ii];
    end
end

if ~isempty(empty_items)
    path_out_tmp = path_out;
    path_out = cell(1, length(path_out_tmp) - length(empty_items));
    cnt = 1;
    cnt1 = 1;
    for ii = 1:length(path_out_tmp)
        if ii ~= empty_items(cnt1)
            path_out{cnt} = path_out_tmp{ii};
            cnt = cnt + 1;
        else
            if cnt1 < length(empty_items)
                cnt1 = cnt1 + 1;
            end
        end
    end
end
