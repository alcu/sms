function [xout yout] = pick_point_clock(x, y)
% function [xout yout] = pick_point_clock(x, y)
%
% Picks the "most clockwise" point out of 2 points. If only 1 point is
% inputted, the same point is outputted. This assumes there is only a pretty
% small angle between points 1 and 2.
%
% Inputs:
%        If 2 input points, each of x and y should be of size [1 2]:
%        - Point 1: x(1) = x-coord, y(1) = y-coord
%        - Point 2: x(2) = x-coord, y(2) = y-coord
%        If 1 intersection point, each of x and y should be of size [1 1]:
%        - x = x-coord, y = y-coord
% Outputs:
%        - xout = x-coord, yout = y-coord

do_plots = 0;
if do_plots
    if length(x) == 2 & length(y) == 2
        figure;
        plot([0 x(1)], [0 y(1)]); hold on
        plot([0 x(2)], [0 y(2)], '-r');
        legend('Point 1', 'Point 2')
    end
end

if length(x) == 1 % Only 1 point inputted.
    xout = x;
    yout = y;
    return
end

ang1 = angle(x(1) + i*y(1));
ang2 = angle(x(2) + i*y(2));

if ang1 > pi/2 & ang2 < -pi/2
    xout = x(1); yout = y(1);
elseif ang2 > pi/2 & ang1 < -pi/2
    xout = x(2); yout = y(2);
else % Points are not around boundary of +/-pi: Pick point with smallest angle.
    if ang1 > ang2
        xout = x(2); yout = y(2);
    else
        xout = x(1); yout = y(1);
    end
end
