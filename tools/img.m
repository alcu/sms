function img(data, varargin)
% function img(data, varargin)
%
% Simple image display function.
%
% Note: This function assumes that your image data matrix, data, is in the
% normal MATLAB format such that the x-position is increasing from up to
% down, and the y-position is increasing from left to right: data(xpos, ypos)
% The function will display the image matrix with a 90-degree
% counterclockwise rotation, such that the x-position is increasing from left
% to right, and the y-position is increasing from down to up.
%
% Inputs:
% data            2D or 3D image matrix. e.g. data(nx, ny, nslcs)
% 'clims'         For imagesc().
% 'cmap'          For colormap().
% 'nrows'         For 3D image matrices only: Number of rows of slices to display.
% 'ncols'         For 3D image matrices only: Number of columns of slices to display.

data = squeeze(data); % Get rid of singleton dimensions.

%%% Defaults.
arg.clims = [];
arg.cmap = gray;
nslcs = size(data, 3);
arg.nrows = floor(sqrt(nslcs));
arg.ncols = ceil(sqrt(nslcs));
if arg.nrows*arg.ncols < nslcs
    arg.nrows = arg.nrows + 1;
end

%%% Replace defaults with user inputs.
arg = vararg_pair(arg, varargin);

[layout order] = tile_rot_layout([arg.nrows arg.ncols]);
data = tile(data, layout, order);

iscomplex = false;
if ~isreal(data)
    iscomplex = true;
    data = abs(data);
    warning('Took abs() of complex data.')
end

if isempty(arg.clims)
    imagesc(data.');
else
    imagesc(data.', arg.clims);
end

axis image;
axis xy; % Flips matrix display up-down. End result (combined with .') is to display the matrix with a 90-degree counterclockwise rotation.
colorbar;
colormap(arg.cmap)

if iscomplex
    title('Abs of complex data')
end
