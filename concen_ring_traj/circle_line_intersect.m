function [x y] = circle_line_intersect(x1, y1, r1, a, b, c)
% function [x y] = circle_line_intersect(x1, y1, r1, a, b, c)
% 
% Finds intersection point(s), if any, of a circle with center (x1, y1)
% and radius r1, and a line with equation y = a*x + b. If the line is vertical, then use [] for a and b, and x = c.
%
% Output:
%        If 2 intersection points, each of x and y will be of size [1 2]:
%        - Point 1: x(1) = x-coord, y(1) = y-coord
%        - Point 2: x(2) = x-coord, y(2) = y-coord
%        If 1 intersection point, each of x and y will be of size [1 1]:
%        - x = x-coord, y = y-coord
%        If no intersection points, x == y == [].

do_plots = 0;
if do_plots
    figure;
    theta = [0:0.01:2*pi];
    plot(x1 + r1*cos(theta), y1 + r1*sin(theta)); hold on
    if isempty(a) & isempty(b)
        ytmp = [-10:0.01:10];
        plot(c*ones(size(ytmp)), ytmp, '-r')
    else
        xtmp = [-10:0.01:10];
        plot(xtmp, a*xtmp + b, '-r')
    end
    axis square
end

% Vertical line.
if isempty(a) & isempty(b)
    x(1) = c;
    b24ac = 4*y1^2 - 4*(c^2 - 2*x1*c + x1^2 + y1^2 - r1^2);
    if b24ac < 0 % No intersection.
        x = []; y = []; return
    else
        y(1) = (2*y1 + sqrt(b24ac))/2;
        if b24ac == 0 % 1 point of intersection.
            return
        else % 2 points of intersection.
            x(2) = c;
            y(2) = (2*y1 - sqrt(b24ac))/2;
            return
        end
    end
else % Non-vertical line.
    b24ac = (-2*x1 + 2*a*b - 2*y1*a)^2 - 4*(1 + a^2)*(x1^2 + b^2 - 2*y1*b + y1^2 - r1^2);
    if b24ac < 0 % No intersection.
        x = []; y = []; return
    else
        x(1) = (2*x1 - 2*a*b + 2*y1*a + sqrt(b24ac))/(2*(1 + a^2));
        y(1) = a*x(1) + b;
        if b24ac == 0 % 1 point of intersection.
            return
        else % 2 points of intersection.
            x(2) = (2*x1 - 2*a*b + 2*y1*a - sqrt(b24ac))/(2*(1 + a^2));
            y(2) = a*x(2) + b;
            return
        end
    end
end
