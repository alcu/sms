function pathstr_out = correct_nfs_link_path(pathstr)
% function pathstr_out = correct_nfs_link_path(pathstr)
%
% If pathstr is a path string to a file/dir, and that file/dir or any
% directory in the full path is a link to somewhere in /net/, and the link is
% broken, the function tries to correct the path and return the corrected
% path in pathstr_out.
%
% Otherwise, it just returns pathstr as pathstr_out.
%
% This function is designed to work even if the file/dir indicated in pathstr
% doesn't exist. This function only works on the first link it finds in
% pathstr starting from the right.
%
% Inputs:
% pathstr: path string containing the file or directory of interest.
%
% Outputs:
% pathstr_out: path string of the corrected location of the file or directory
% of interest.

path_names = cut_path(pathstr);

if pathstr(1) == filesep
    beg_filesep = 1;
    path_clean = fullfile(filesep, fullfile_cell(path_names));
else
    beg_filesep = 0;
    path_clean = fullfile_cell(path_names);
end
if pathstr(end) == filesep
    end_filesep = 1;
else
    end_filesep = 0;
end

if exist(path_clean, 'file') % Don't do anything if the path to the file/dir exists.
    return_orig = 1;
else
    return_orig = 0;
    name = cell(1, length(path_names));
    [status cmdout] = system(['readlink ' escape_fn(path_clean)]);
    if isempty(cmdout) % pathstr does not point to a link.
        if length(path_names) > 1
            for ii = [length(path_names) : -1 : 2]
                name{ii} = path_names{ii};
                path_names{ii} = []; % Check parent file part.
                if beg_filesep
                    path_clean = fullfile(filesep, fullfile_cell(path_names));
                else
                    path_clean = fullfile_cell(path_names);
                end
                [status cmdout] = system(['readlink ' escape_fn(path_clean)]);
                if ~isempty(cmdout) % pathstr points to a link.
                    if ~exist(path_clean, 'file') % Keep on going if the path to the file/dir exists.
                        break
                    end
                elseif ii == 2
                    return_orig = 1;
                end
            end
        else
            return_orig = 1;
        end
    end
end

if return_orig
    pathstr_out = pathstr;
else
    link_path = cmdout(1:end - 1); % Last character of cmdout will be a carriage return.
    link_parts = cut_path(link_path);
    if isequal(link_parts{1}, 'net')
        if isequal(link_parts{3}, 'export')
            link_parts{3} = [];
            pathstr_out = fullfile(filesep, fullfile_cell(link_parts), fullfile_cell(name));
        else
            hname = link_parts{2};
            link_parts{1} = [];
            link_parts{2} = [];
            pathstr_out = fullfile(filesep, 'net', hname, 'export', fullfile_cell(link_parts), fullfile_cell(name));
        end
        if end_filesep
            pathstr_out = [pathstr_out filesep];
        end
    else
        pathstr_out = pathstr;
    end
end
