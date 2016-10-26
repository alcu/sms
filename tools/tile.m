function tmap = tile(varargin)
% function tmap = tile(maps, optional_layout, optional_order)
% Produces a tiled 2D matrix from a 3D matrix.
% Inputs:
%     Required:
%         3D matrix: maps(FOVx, FOVy, slice_num)
%     Optional:
%         optional_layout: tile layout -> [num_rows num_columns]
%             If this optional input is left out, this function will tile the
%             3D matrix in the smallest possible square 2D matrix.
%         optional_order: e.g. [4 1 5 2 6 3] for a 3-by-2 layout.
%             If this optional input is left out, the order is just [1:num_slices].
% 
% Output: tiled 2D matrix.

maps = varargin{1};

xsize = size(maps, 1);
ysize = size(maps, 2);
zsize = size(maps, 3);

if nargin == 1
    sqsize = ceil(sqrt(zsize));
    
    repnum = floor(zsize/sqsize);
    
    extranum = zsize - repnum*sqsize;
    
    temp = zeros(sqsize*xsize, sqsize*ysize);
    
    for rownum = 1:repnum
        %size(temp(((1:xsize) + (rownum - 1)*xsize), :))
        %size(maps(:, :, ((1:sqsize) + ((rownum - 1)*sqsize))))
        
        temp(((1:xsize) + (rownum - 1)*xsize), :) = reshape(maps(:, :, ((1:sqsize) + ((rownum - 1)*sqsize))), [xsize ysize*sqsize]);
    end;
    %disp('hi')
    for xnum = 1:extranum
        %disp('hi')
        temp(((1:xsize) + repnum*xsize), ((1:ysize) + ((xnum - 1)*ysize))) = maps(:, :, ((repnum*sqsize) + xnum));
    end;
    tmap = temp;
else
    layout = varargin{2};
    if prod(layout) < zsize
        error('Tiling layout too small.')
    end
    if nargin == 2
        ord = [1:zsize];
    else
        ord = varargin{3};
    end
    if length(ord) < zsize
        error('Length of order vector is smaller than the number of slices.')
    end
    
    tmap = zeros(layout(1)*xsize, layout(2)*ysize);
    slc_cnt = 0;
    max_cnt = zsize;
    for rw = 1:layout(1)
        for cl = 1:layout(2)
            slc_cnt = slc_cnt + 1;
            if ord(slc_cnt) > zsize
                max_cnt = max_cnt + 1;
                continue;
            end
            tmap((rw - 1)*xsize + 1 : rw*xsize, ...
                 (cl - 1)*ysize + 1 : cl*ysize) = maps(:, :, ord(slc_cnt));
            if slc_cnt == max_cnt
                break;
            end
        end
        if slc_cnt == max_cnt
            break;
        end
    end
end
