function [layoutfl order] = tile_rot_layout(layout)
% function [layoutfl order] = tile(layout)
% Produces desired layout and order for use with tile() if you want to rotate
% the matrix so that the x-position is increasing from left to right, and the
% y-position is increasing from down to up.
%
% Inputs:
%     Required:
%         layout -> [num_rows num_columns]
% 
% Outputs: 
%     layoutfl and order. Use with tile().

layoutfl = layout(2:-1:1);

order = reshape([1:prod(layout)], layoutfl);
order = flipud(order.');

order = order(:).';
